
# TF-Gene Regulatory Network Pipeline

![Python](https://img.shields.io/badge/python-3.8%2B-blue)

This repository provides a **notebook-aligned Python pipeline** to construct transcription factor (TF)–gene regulatory networks from tumor and normal expression data. It performs DEG (Differentially Expressed Genes) analysis, correlation-based network inference, Ridge regression-based network modeling, network cleaning, visualization, and summary generation.

---

##  Features

* Combines tumor and normal expression datasets into a unified matrix
* Filters TFs and computes DEGs
* Infers:

  * Correlation-based TF-gene network
  * Ridge regression-based TF-gene network (supports GPU via cuML)
* Cleans network using top-k and threshold filters
* Generates static plots and interactive HTML visualizations
* Provides summary of master regulators and strongest interactions
* Notebook-aligned workflow for reproducibility

---

##  Repository Structure

```
## Project Structure

tf_gene_network/
│
├── data/                       # Only smaller or processed files
│   ├── tf_list.csv             # Optional TF list derived from TXT
│   ├── combined_expression.csv  # Combined expression matrix (auto-generated)
│   ├── DEG_filtered.csv        # Differentially expressed genes (auto-generated)
│   └── ...                     # Other intermediate files
│
├── results/                    # Output directory for plots and network files
│   ├── correlation_distribution.png
│   ├── top_tfs_bar.png
│   ├── correlation_small_network.png
│   ├── clean_tf_gene_network.png
│   ├── interactive_tf_gene_network.html
│   ├── tf_gene_network_raw.csv
│   ├── tf_gene_network_clean.csv
│   └── network_summary.txt
│
├── src/                        # Source code
│   └── tf_gene_network.py       # Main pipeline script
│
├── notebooks/                   # Jupyter notebooks for analysis
│   └── ...                     
├── requirements.txt             # Python dependencies
├── .gitignore                   # Files/folders to ignore in Git
└── README.md                    # Project documentation
```
### Datasets

> Large raw expression files are not included to reduce repo size. Download and place them in `data/`.

| File | Description | Link |
|------|-------------|------|
| `BC_TCGA-Normal.txt` | Normal breast tissue expression | [Download](https://drive.google.com/drive/folders/1cupQ78ABMrqFRLEhm9vua8_FUoqFA6iT) |
| `BC_TCGA-Tumor.txt` | Tumor breast tissue expression | [Download](https://drive.google.com/drive/folders/1cupQ78ABMrqFRLEhm9vua8_FUoqFA6iT) |
| `combined_expression.csv` | Tumor + Normal expression | [Download](https://drive.google.com/drive/folders/1cupQ78ABMrqFRLEhm9vua8_FUoqFA6iT) |

##  Installation```
```
```
1. **Clone the repository**

```bash
git clone repo_url.git
cd tf_gene_network
```

2. **Create a virtual environment**

```bash
python -m venv venv
source venv/bin/activate       # Linux/macOS
venv\Scripts\activate          # Windows
```

3. **Install dependencies**

```bash
pip install -r requirements.txt
```

**Optional GPU support**:

```bash
pip install cupy-cuda12x cuml-cuda12x
```

4. **Dependencies**

* numpy
* pandas
* matplotlib
* networkx
* scipy
* statsmodels
* scikit-learn
* pyvis (for interactive HTML visualization)
* cupy + cuml (optional for GPU Ridge regression)

---

##  Usage

Run the entire pipeline:

```bash
python tf_gene_network.py
```

**Default directories**:

* Input data: `data/`
* Output results: `results/`

**Pipeline steps**:

1. Load or create combined expression matrix
2. Load TF list
3. Filter TF expression
4. Load or compute DEGs
5. Prepare TF and target sets
6. Correlation baseline network
7. Ridge regression network
8. Clean network
9. Save network files
10. Generate plots (correlation distribution, top TFs, small network)
11. Plot cleaned network
12. Create interactive HTML network (if `pyvis` installed)
13. Summarize network

---

##  Outputs

1. **Combined expression**: `combined_expression.csv`
2. **Filtered TF expression**: `filtered_TF_expression.csv`
3. **DEGs**: `DEG_filtered.csv`
4. **Networks**:

   * Raw Ridge network: `tf_gene_network_raw.csv`
   * Cleaned Ridge network: `tf_gene_network_clean.csv`
   * Correlation network: `TF_network_correlation.csv`
5. **Plots**:

   * Correlation distribution: `correlation_distribution.png`
   * Top TF bar chart: `top_tfs_bar.png`
   * Small network: `small_network.png`
   * Clean TF-gene network: `clean_tf_gene_network.png`
6. **Interactive HTML**: `interactive_tf_gene_network.html`
7. **Network summary**: `network_summary.txt`

---

##  Configuration

The script allows configuration through constants and function parameters:

```python
BASE_PATH = "data"           # Input directory
OUT_DIR = "results"          # Output directory
CORR_THRESHOLD = 0.7         # Correlation network threshold
RIDGE_COEF_THRESH = 1e-3     # Ridge network coefficient threshold
CLEAN_K = 3                  # Top-k edges per target
CLEAN_THRESHOLD = 0.05       # Minimum weight for cleaned network
```

Modify `run_pipeline()` arguments to change behavior.

---
### Results / Output Files

> Some large network CSVs (`tf_gene_network_ridge.csv`, `tf_gene_network_raw.csv`) are **not included** due to size. The pipeline generates them in `results/` when you run it.

| File | Description |
|------|-------------|
| `correlation_distribution.png` | Correlation histogram of TF vs target genes |
| `top_tfs_bar.png` | Bar chart of top TFs (correlation) |
| `correlation_small_network.png` | Small network plot (correlation) |
| `clean_tf_gene_network.png` | Cleaned TF-gene network plot |
| `interactive_tf_gene_network.html` | Interactive network visualization |
| `network_summary.txt` | Summary of raw & clean networks |
| `tf_gene_network_ridge.csv` | Raw Ridge network (large) – **generated by pipeline** |
| `tf_gene_network_raw.csv` | Cleaned Ridge network (large) – **generated by pipeline** |
##  Notes

* Ensure your expression files have **genes as rows** and **samples as columns**
* TF list file can be `Homo_sapiens_TF.txt` or `tf_list.csv` (preferred)
* GPU Ridge regression is optional; CPU fallback is automatically used
* PyVis is required only for interactive HTML network visualization


