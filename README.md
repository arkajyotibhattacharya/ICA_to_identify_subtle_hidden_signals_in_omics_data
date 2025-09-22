
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


---

## References

- Bhattacharya, A., Bense, R.D., Urzúa-Traslaviña, C.G. **et al.** (2020). *Transcriptional effects of copy number alterations in a large set of human cancers.* **Nature Communications**, 11, 715. https://doi.org/10.1038/s41467-020-14605-5

- Knapen, D.G., Hone Lopez, S., de Groot, D.J.A. **et al.** (2024). *Independent transcriptional patterns reveal biological processes associated with disease-free survival in early colorectal cancer.* **Communications Medicine**, 4, 79. https://doi.org/10.1038/s43856-024-00504-z

- Hone Lopez, S., Loipfinger, S., Bhattacharya, A. **et al.** (2025). *Upfront whole blood transcriptional patterns in patients receiving immune checkpoint inhibitors associate with clinical outcome.* **Cancer Immunology, Immunotherapy**, 74, 301. https://doi.org/10.1007/s00262-025-04155-4

- Leeuwenburgh, V.C., Urzúa-Traslaviña, C.G., Bhattacharya, A. **et al.** (2021). *Robust metabolic transcriptional components in 34,494 patient-derived cancer-related samples and cell lines.* **Cancer & Metabolism**, 9, 35. https://doi.org/10.1186/s40170-021-00272-7

- Bhattacharya, A., Stutvoet, T.S., Perla, M. **et al.** (2025). *Transcriptional pattern enriched for synaptic signaling is associated with shorter survival of patients with high-grade serous ovarian cancer.* **eLife** (Reviewed Preprint), Version of Record May 13, 2025. https://doi.org/10.7554/eLife.101369.2

- Urzúa-Traslaviña, C.G., Leeuwenburgh, V.C., Bhattacharya, A. **et al.** (2021). *Improving gene function predictions using independent transcriptional components.* **Nature Communications**, 12, 1464. https://doi.org/10.1038/s41467-021-21671-w

