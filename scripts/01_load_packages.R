# Script: 01_load_packages.R
# Project: PCOS WGCNA (Biomedicines 2023)
# Author: Dr. Roozbeh Heidarzadehpilehrood
#
# Minimal helper to load (and, if needed, install) the packages I use
# across the PCOS WGCNA analysis. Keep this file short and boring.

packages <- c(
  "GEOquery",
  "affy",
  "limma",
  "WGCNA",
  "dynamicTreeCut",
  "tidyverse",
  "Biobase"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

# End of 01_load_packages.R
