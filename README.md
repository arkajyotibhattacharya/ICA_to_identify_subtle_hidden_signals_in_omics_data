# ICA Expression Analysis Pipeline (Shiny App)

This repository contains a Shiny app that performs **Independent Component Analysis (ICA)** on gene expression data, with optional clinical annotation and survival analysis.

---

## 🚀 Pipeline Overview

The pipeline runs:

1. **QC** → summary stats + histograms  
2. **PCA** → checks for background effect (removes PC1 if needed)  
3. **ICA** → decomposes expression data into independent components  
4. **IC exploration** → scatter plots, optional GSEA  
5. **Associations** → with sample annotations  
6. **Survival tree** → if survival time + status are provided  

---

## 📥 Input files

- **Expression profile (required)**: genes × samples table (TSV/CSV).  
- **Genomic mapping file (optional)**: maps probes → Entrez/symbol + chr/bp.  
- **Sample annotation file (optional)**: sample metadata (rows = samples).  
- **Survival columns (optional, inside annotation)**:  
  - Survival time (e.g. `OS`)  
  - Survival status (e.g. `OS.binary`, 0 = censored, 1 = event)  
- **Gene set GMT file (optional)**: for GSEA.

---

## ⚙️ What the app does

- Deduplicates IDs (randomly keeps one if duplicated).  
- Converts probes/symbols → Entrez IDs if mapping is provided.  
- Sorts genes by chromosome and base pair if mapping is provided.  
- Skips steps automatically if files are missing (e.g. no annotation → no association analysis).  

---

## 📤 Outputs

All results are saved under `Results/` and `Data/`.  

- `expr_qc_genes_hist.pdf` / `expr_qc_samples_hist.pdf`  
- `pca_summary.tsv` (variance explained + PC1 info)  
- `independent_components.txt` (genes × ICs)  
- `mixing_matrix.txt` (samples × ICs)  
- `independent_components_scatter.pdf`  
- `gsea_heatmap.pdf` (if GSEA run)  
- `associations_heatmap.pdf` (if annotations provided)  
- `survival_tree.pdf` (if survival columns provided)

---

## ▶️ Running the app

Clone this repo and run in R:

```r
install.packages(c("shiny","data.table","fastICA","survival","ggplot2"))
# optional:
# install.packages(c("pheatmap","MST"))

shiny::runApp("app.R")

# ICA_to_identify_subtle_hidden_signals_in_omics_data
