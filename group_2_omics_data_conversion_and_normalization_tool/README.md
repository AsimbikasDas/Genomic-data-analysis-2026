# 🧬 OmicsForge — Ultra Fast RNA Normalisation

**OmicsForge** is a State-Of-The-Art (SOTA) bioinformatics web application designed to instantly convert raw RNA-Sequencing read counts into mathematically normalized metrics: **TPM (Transcripts Per Million)** and **RPKM (Reads Per Kilobase Million)**. 

By automating the correction for both *Sequencing Depth Bias* and *Gene Length Bias*, this tool enables researchers to perform fair, publication-ready gene expression comparisons in milliseconds.

---

## 🏗️ System Architecture

The project has been decoupled into a modern MVC architecture for maximum performance and a beautiful user experience.

*   **Frontend (`/frontend`)**: A gorgeous, highly responsive React application built with **Next.js** and **Tailwind CSS**. It utilizes a Glassmorphism design system, abstract DNA-themed background artwork, and Framer Motion micro-interactions.
*   **Backend (`/backend`)**: A lightning-fast RESTful API powered by **FastAPI**. It handles complex matrix operations using **Pandas** and **NumPy** vectorization, ensuring extremely large CSV files are processed almost instantaneously without blocking the main event loop.

---

## 🚀 How to Run the Project

To launch the project locally, you will need to start both the Python Backend and the Next.js Frontend in two separate terminal windows.

### 1. Start the Backend API (Terminal 1)
```powershell
# Navigate to the backend directory
cd backend

# Install the required Python dependencies
pip install -r requirements.txt

# Launch the FastAPI Server
python -m uvicorn main:app --reload --port 8000
```
*The backend will now be actively listening at `http://localhost:8000`.*

### 2. Start the Frontend UI (Terminal 2)
```powershell
# Navigate to the frontend directory
cd frontend

# Install the strict Node dependencies (only required the first time)
npm install

# Launch the Next.js Development Server
npm run dev
```
*The frontend will launch at `http://localhost:3000` (or `3001` if port 3000 is occupied). Simply visit this link in your web browser to use OmicsForge.*

---

## 👥 Project Members (Group 2)

| Project Member Name | Roll Number |
| :--- | :--- |
| **S Kedareswar** | 22BTB0A37 |
| **Suryaansh Dev** | 22BTB0A76 |
| **Lokesh Verma** | 22BTB0A63 |
| **Anugu Nithish** | 22BTB0A32 |
