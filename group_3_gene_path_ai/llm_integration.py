import os
from huggingface_hub import InferenceClient
from dotenv import load_dotenv 

load_dotenv() 

import os
from huggingface_hub import InferenceClient

def generate_biological_summary(up_regulated, down_regulated, pathways_df, tissue_type="Unknown", detail_level="Concise Summary"):
    # 1. Initialize the client with Llama-3
    client = InferenceClient(
        model="meta-llama/Meta-Llama-3-8B-Instruct", 
        token=os.environ.get("HF_TOKEN") 
    )
    
    # 2. Safely extract data (prevents double-joining and JSON breaking)
    if not isinstance(pathways_df, str) and not pathways_df.empty:
        top_pathways = [str(term).replace('"', "'") for term in pathways_df['Term'].head(10).tolist()]
    else:
        top_pathways = ["No pathways found"]
        
    safe_up = [str(g).replace('"', "'") for g in up_regulated[:30]]
    safe_down = [str(g).replace('"', "'") for g in down_regulated[:30]]
    
    # 3. Determine specific instructions based on detail level
    if detail_level == "Concise Summary":
        instruction = "Provide a highly concise summary focusing on the immediate biological significance."
    elif detail_level == "Detailed Mechanistic Report":
        instruction = "Provide a detailed step-by-step explanation of the molecular mechanisms and how these pathways interact."
    else: 
        instruction = "Focus heavily on the translational impact. Explain how these pathway shifts could be targeted therapeutically or used as clinical biomarkers."

    # 4. Construct the highly-structured prompt
    prompt = f"""
    You are an expert Computational Biologist analyzing RNA-seq differential expression data.
    
    DATA CONTEXT:
    - Experimental Condition / Tissue: {tissue_type}
    - Top Up-regulated Genes: {', '.join(safe_up)}
    - Top Down-regulated Genes: {', '.join(safe_down)}
    - Top Enriched Pathways: {', '.join(top_pathways)}
    
    YOUR TASK:
    Write a cohesive biological interpretation of what is fundamentally happening in these cells. 
    {instruction}
    
    STRICT RULES:
    1. DO NOT just list the pathways or repeat the data back to me. 
    2. DO NOT write "Pathway X (genes A, B, C)". 
    3. Synthesize the data into "Biological Themes".
    4. Explain the *mechanism*. Keep in mind that Up-regulated means higher in the Experimental group, and Down-regulated means higher in the Baseline group.
    5. Avoid assuming cancer unless the data or tissue type explicitly suggests it.
    
    FORMAT YOUR RESPONSE WITH THESE HEADINGS:
    ### 🔬 Executive Biological Summary
    ### ⚙️ Key Mechanistic Drivers
    ### 🩺 Phenotypic & Clinical Implications
    """
    
    try:
        # BUG FIX: Added "role": "user" to fix the Bad Request error
        messages = [{"role": "user", "content": prompt}]
        
        # 5. Call the chat completion with Llama-3
        response = client.chat_completion(
            model="meta-llama/Meta-Llama-3-8B-Instruct", 
            messages=messages,
            max_tokens=1000, 
            temperature=0.2
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Error connecting to LLM: {e}"