# Script: 02_download_and_preprocess_GSE48301.R
# Project: PCOS WGCNA (Biomedicines 2023)
# Author: Dr. Roozbeh Heidarzadehpilehrood
#
# Goal:
#   - Download GSE48301 from GEO
#   - Build expression matrix and phenotype table
#   - Apply a light expression filter (remove completely dead probes)
#   - Save processed objects for downstream limma + WGCNA steps

source("scripts/01_load_packages.R")
options(stringsAsFactors = FALSE)

gse_id <- "GSE48301"

# Output files live in data/
if (!dir.exists("data")) dir.create("data")

expr_outfile_raw  <- file.path("data", paste0(gse_id, "_expr_raw.rds"))
expr_outfile_filt <- file.path("data", paste0(gse_id, "_expr_filtered.rds"))
pheno_outfile     <- file.path("data", paste0(gse_id, "_pheno.csv"))

# -------------------------------------------------------------------
# 1) Download series and pull out the ExpressionSet
# -------------------------------------------------------------------

gset <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
if (length(gset) > 1) {
  eset <- gset[[1]]  # GSE48301 only has one anyway, but just in case
} else {
  eset <- gset[[1]]
}

expr_mat <- Biobase::exprs(eset)
pheno_df <- Biobase::pData(eset)

cat("Raw expression matrix dimensions (probes x samples):\n")
print(dim(expr_mat))
cat("Raw phenotype dimensions:\n")
print(dim(pheno_df))

# -------------------------------------------------------------------
# 2) Basic sanity check on expression scale
# -------------------------------------------------------------------
# Note: this dataset is already log-ish; I am not trying to be fancy here.
# If I ever re-use this script on another array, I can revisit this part.

summary(as.numeric(expr_mat))

# -------------------------------------------------------------------
# 3) Light expression filtering
# -------------------------------------------------------------------
# Idea: drop probes that are essentially off in almost all samples.
# I do NOT want to end up with only a few dozen genes for WGCNA.

expr_threshold <- 3      # intensity cutoff
min_prop       <- 0.20   # must be "on" in at least 20% of samples

keep_probes <- rowMeans(expr_mat > expr_threshold, na.rm = TRUE) >= min_prop
table(keep_probes)

expr_filt <- expr_mat[keep_probes, , drop = FALSE]

cat("Genes before filtering:", nrow(expr_mat),  "\n")
cat("Genes after  filtering:", nrow(expr_filt), "\n")

# If this ever collapses to a tiny number again, I can temporarily skip filtering:
# expr_filt <- expr_mat

# -------------------------------------------------------------------
# 4) Save processed objects
# -------------------------------------------------------------------

saveRDS(expr_mat,  expr_outfile_raw)
saveRDS(expr_filt, expr_outfile_filt)
write.csv(pheno_df, pheno_outfile, row.names = TRUE)

message("Saved raw expression matrix to:  ", expr_outfile_raw)
message("Saved filtered expression matrix to: ", expr_outfile_filt)
message("Saved phenotype data to: ", pheno_outfile)

# End of 02_download_and_preprocess_GSE48301.R
