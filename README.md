README.md

# PCOS WGCNA – Biomedicines 2023

This repository contains the reproducible R code used in the manuscript:

> Heidarzadehpilehrood R, Pirhoushiaran M, Osman MB, Hamid HA, Ling KH.  
> **Weighted gene co-expression network analysis (WGCNA) discovered novel long non-coding RNAs for polycystic ovary syndrome.**  
> *Biomedicines* 2023;11(2):518.

The goal of this project is to:
- Download and preprocess the public microarray dataset **GSE48301** (PCOS vs controls)
- Perform differential expression analysis (PCOS vs control)
- Build a WGCNA network and detect co-expression modules
- Correlate modules with the PCOS trait
- Export gene–module assignments, module–trait correlations, and example figures

All analyses are implemented in R and organised as a small, self-contained pipeline.

---

## Directory structure

- `pcos-wgcna-biomedicines-2023.Rproj`  
  RStudio project file (recommended entry point).

- `scripts/`  
  Main analysis scripts:
  - `01_load_packages.R`  
    Install and load required packages (GEOquery, affy, limma, WGCNA, tidyverse, etc.).
  - `02_download_and_preprocess_GSE48301.R`  
    Download GSE48301 from GEO, extract the expression matrix and phenotype table, apply basic filtering, and save:
    - `data/GSE48301_expr_raw.rds`
    - `data/GSE48301_expr_filtered.rds`
    - `data/GSE48301_pheno.csv`
  - `03_DE_analysis_limma_DEGs_DElncR.R`  
    Limma differential expression analysis (PCOS vs control) and export of the DEG table to:
    - `results/GSE48301_DEG_PCOS_vs_Control_limma.csv`
  - `04_WGCNA_network_and_modules_GSE48301.R`  
    WGCNA network construction, module detection, module–trait correlation, and export of:
    - `results/GSE48301_WGCNA_gene_module_assignment.csv`
    - `results/GSE48301_WGCNA_module_trait_correlations.csv`  
    as well as R objects in `data/` for downstream analyses.

- `data/`  
  Processed expression matrices, phenotype table, and saved WGCNA objects.

- `data-raw/`  
  Optionally used for original GEO downloads (series matrix text file), if needed.

- `results/`  
  CSV outputs from the DE and WGCNA analyses.

- `figures/`  
  (Optional) Plots generated from the pipeline (e.g. sample clustering, scale-free topology, module–trait heatmap).

---

## Requirements

- R (≥ 4.2 recommended)
- R packages:
  - `GEOquery`
  - `affy`
  - `limma`
  - `WGCNA`
  - `Biobase`
  - `tidyverse`
  - `dynamicTreeCut`

All required packages are loaded (and installed if missing) by `scripts/01_load_packages.R`.

---

## How to run the pipeline

1. Clone this repository or download it as a ZIP and extract it.
2. Open `pcos-wgcna-biomedicines-2023.Rproj` in RStudio.
3. Run the scripts in the following order, linearly:
   1. `scripts/01_load_packages.R`
   2. `scripts/02_download_and_preprocess_GSE48301.R`
   3. `scripts/03_DE_analysis_limma_DEGs_DElncR.R`
   4. `scripts/04_WGCNA_network_and_modules_GSE48301.R`
4. Outputs will be written to the `data/` and `results/` directories.  
   Figures can be generated from the WGCNA objects and saved to `figures/`.

---

   ## License
Code in this repository is released under the MIT License.  
The underlying expression data (GSE48301) remain subject to the original GEO / authors’ terms of use.

---

## Author

**Dr Roozbeh Heidarzadehpilehrood**  
Human Geneticist – Transcriptomics & ncRNA biomarkers  
Contact: heidarzadeh.roozbeh [at] gmail [dot] com
