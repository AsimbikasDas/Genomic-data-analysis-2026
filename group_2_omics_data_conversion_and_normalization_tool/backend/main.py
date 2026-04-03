import os
import io
import pandas as pd
from fastapi import FastAPI, UploadFile, File, Form, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse
from fastapi.exceptions import RequestValidationError
from starlette.exceptions import HTTPException as StarletteHTTPException

from services import parse_csv_file, get_preview_data, process_normalization
from analysis import run_full_analysis

app = FastAPI(title="OmicsForge API", description="SOTA RNA-Seq Normalization Backend")

# Dynamic CORS configuration
raw_origins = os.environ.get("ALLOWED_ORIGINS", "")
allowed_origins = [o.strip().rstrip("/") for o in raw_origins.split(",") if o.strip()]

# Hardcoded fallbacks for deployment
for origin in ["https://project-roan-six-31.vercel.app", "http://localhost:3000"]:
    if origin not in allowed_origins:
        allowed_origins.append(origin)

# If no restricted environment variable is set, allow all for debugging/testing
if not raw_origins:
    allowed_origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=True if "*" not in allowed_origins else False,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    return JSONResponse(status_code=400, content={"detail": f"Internal Error: {str(exc)}"})

@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    return JSONResponse(status_code=422, content={"detail": f"Validation Error: {str(exc)}"})

@app.exception_handler(StarletteHTTPException)
async def http_exception_handler(request: Request, exc: StarletteHTTPException):
    return JSONResponse(status_code=exc.status_code, content={"detail": exc.detail})

@app.get("/")
def health_check():
    return {"status": "OmicsForge API is Online", "cors_mode": "Permissive (Wildcard)" if "*" in allowed_origins else "Restricted"}

@app.post("/api/preview")
async def preview_csv(file: UploadFile = File(...)):
    """Receive file payloads, dispatch to Services, and return metadata properties."""
    if not file.filename.endswith(".csv"):
        raise HTTPException(status_code=400, detail="Invalid file format! Please upload a CSV.")
    
    try:
        content = await file.read()
        df = parse_csv_file(content)
        result = get_preview_data(df)
        return JSONResponse(content=result)
    except Exception as e:
        return JSONResponse(status_code=400, content={"detail": str(e)})

@app.post("/api/normalize")
async def normalize_csv(
    file: UploadFile = File(...),
    gene_id_col: str = Form(...),
    do_tpm: str = Form("true"),
    do_rpkm: str = Form("true")
):
    """Receive files, execute computational transformations, and stream normalized datasets back as CSV."""
    try:
        is_tpm = do_tpm.lower() == "true"
        is_rpkm = do_rpkm.lower() == "true"
        
        content = await file.read()
        df = parse_csv_file(content)
        
        # Now returns a DataFrame thanks to our services.py update
        normalized_df = process_normalization(df, gene_id_col, is_tpm, is_rpkm)
        
        # Truly stream the CSV in chunks to save memory for very large datasets
        def csv_generator():
            # first chunk including header
            chunk_size = 1000
            for i in range(0, len(normalized_df), chunk_size):
                # Only include header for the first chunk
                yield normalized_df.iloc[i : i + chunk_size].to_csv(index=False, header=(i == 0))
            
        return StreamingResponse(
            csv_generator(), 
            media_type="text/csv",
            headers={"Content-Disposition": "attachment; filename=normalized_data.csv"}
        )
        
    except Exception as e:
        print(f"Normalisation Error: {str(e)}")
        return JSONResponse(status_code=400, content={"detail": str(e)})

@app.post("/api/analyze")
async def analyze_data(
    file: UploadFile = File(...),
    gene_id_col: str = Form(...),
    method: str = Form("spearman")
):
    """Run analysis pipelines (correlation, PCA, distributions) on normalized data."""
    try:
        content = await file.read()
        df = parse_csv_file(content)
        results = run_full_analysis(df, gene_id_col, method)
        return JSONResponse(content=results)
    except Exception as e:
        print(f"Analysis Error: {str(e)}")
        return JSONResponse(status_code=400, content={"detail": str(e)})
