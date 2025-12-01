# Script: 03_DE_analysis_limma_DEGs_DElnc.R
# Project: PCOS WGCNA (Biomedicines 2023)
# Author: Dr. Roozbeh Heidarzadehpilehrood
#
# Goal:
#   - Load filtered expression + phenotype data for GSE48301
#   - Build PCOS vs Control contrast
#   - Run limma and export a DEG table

source("scripts/01_load_packages.R")
options(stringsAsFactors = FALSE)

gse_id <- "GSE48301"

expr_outfile_filt <- file.path("data",   paste0(gse_id, "_expr_filtered.rds"))
pheno_outfile     <- file.path("data",   paste0(gse_id, "_pheno.csv"))
deg_outfile       <- file.path("results", paste0(gse_id, "_DEG_PCOS_vs_Control_limma.csv"))

if (!dir.exists("results")) dir.create("results")

# -------------------------------------------------------------------
# 1) Load data
# -------------------------------------------------------------------

expr_filt <- readRDS(expr_outfile_filt)
pheno_df  <- read.csv(pheno_outfile, row.names = 1, check.names = FALSE)

cat("Filtered expression matrix dimensions (genes x samples):\n")
print(dim(expr_filt))
cat("Phenotype dimensions:\n")
print(dim(pheno_df))

# -------------------------------------------------------------------
# 2) Match samples between expression and phenotype
# -------------------------------------------------------------------

common_samples <- intersect(colnames(expr_filt), rownames(pheno_df))
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between expression and phenotype.")
}

expr_filt <- expr_filt[, common_samples]
pheno_df  <- pheno_df[common_samples, ]

cat("After matching samples:\n")
print(dim(expr_filt))
print(dim(pheno_df))

# -------------------------------------------------------------------
# 3) Define PCOS vs Control groups
# -------------------------------------------------------------------
# Strategy: find a phenotype column that actually contains the string "PCOS".

is_text_col    <- sapply(pheno_df, function(x) is.character(x) || is.factor(x))
candidate_cols <- names(pheno_df)[is_text_col]

pcos_cols <- names(pheno_df)[sapply(
  pheno_df[, candidate_cols, drop = FALSE],
  function(col) any(grepl("PCOS", col, ignore.case = TRUE))
)]

if (length(pcos_cols) == 0) {
  stop("Could not find any phenotype column containing 'PCOS'.")
}

diag_col <- pcos_cols[1]
cat("Using diagnosis column for PCOS vs Control:", diag_col, "\n")
cat("Example values:\n")
print(head(pheno_df[[diag_col]]))

group <- ifelse(
  grepl("PCOS", pheno_df[[diag_col]], ignore.case = TRUE),
  "PCOS",
  "Control"
)
group <- factor(group, levels = c("Control", "PCOS"))

cat("Group counts:\n")
print(table(group))

# -------------------------------------------------------------------
# 4) limma model: PCOS vs Control
# -------------------------------------------------------------------

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "PCOS_vs_Control")

fit  <- limma::lmFit(expr_filt, design)
fit2 <- limma::eBayes(fit)

deg_table <- limma::topTable(
  fit2,
  coef  = "PCOS_vs_Control",
  number = Inf,
  sort.by = "P"
)

head(deg_table)

# -------------------------------------------------------------------
# 5) Export DEG table
# -------------------------------------------------------------------

write.csv(deg_table, deg_outfile, row.names = TRUE)
message("Saved DEG table to: ", deg_outfile)

# End of 03_DE_analysis_limma_DEGs_DElnc.R
