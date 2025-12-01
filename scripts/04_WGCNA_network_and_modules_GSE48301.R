# Script: 04_WGCNA_network_and_modules_GSE48301.R
# Project: PCOS WGCNA (Biomedicines 2023)
# Author: Dr. Roozbeh Heidarzadehpilehrood
#
# Goal:
#   - Take the filtered GSE48301 expression matrix
#   - Construct a WGCNA network
#   - Identify modules and correlate them with PCOS vs Control
#   - Save gene–module assignment, module–trait correlations,
#     and a few key figures as PNGs for the GitHub repo

source("scripts/01_load_packages.R")
options(stringsAsFactors = FALSE)

library(Biobase)
library(dynamicTreeCut)
WGCNA::allowWGCNAThreads()

gse_id <- "GSE48301"

expr_outfile_filt <- file.path("data",   paste0(gse_id, "_expr_filtered.rds"))
pheno_outfile     <- file.path("data",   paste0(gse_id, "_pheno.csv"))

if (!dir.exists("results"))  dir.create("results")
if (!dir.exists("figures"))  dir.create("figures")

# -------------------------------------------------------------------
# 1) Load processed data
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

expr_filt <- expr_filt[, common_samples, drop = FALSE]
pheno_df  <- pheno_df[common_samples, , drop = FALSE]

cat("After matching samples:\n")
print(dim(expr_filt))
print(dim(pheno_df))

# -------------------------------------------------------------------
# 3) Build PCOS vs Control trait (same logic as in script 03)
# -------------------------------------------------------------------

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

traitPCOS <- ifelse(group == "PCOS", 1, 0)
traitData <- data.frame(PCOSvsControl = traitPCOS)
rownames(traitData) <- common_samples

# -------------------------------------------------------------------
# 4) Select top variable genes for WGCNA
# -------------------------------------------------------------------
# Keep up to 10k most variable genes: enough to see structure,
# still manageable on a laptop.

n_top_genes <- min(10000, nrow(expr_filt))

gene_var   <- apply(expr_filt, 1, var, na.rm = TRUE)
ord        <- order(gene_var, decreasing = TRUE)
keep       <- ord[1:n_top_genes]
expr_wgcna <- expr_filt[keep, , drop = FALSE]

cat("Keeping top", n_top_genes, "most variable genes for WGCNA.\n")
cat("Expression matrix for WGCNA (genes x samples):\n")
print(dim(expr_wgcna))

# WGCNA wants samples in rows and genes in columns
datExpr0 <- t(expr_wgcna)
cat("datExpr0 dimensions (samples x genes):\n")
print(dim(datExpr0))

# -------------------------------------------------------------------
# 5) goodSamplesGenes check
# -------------------------------------------------------------------

gsg <- WGCNA::goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    cat("Removing bad genes:",
        paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", "), "\n")
  if (sum(!gsg$goodSamples) > 0)
    cat("Removing bad samples:",
        paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "), "\n")
  
  datExpr0  <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  traitData <- traitData[rownames(datExpr0), , drop = FALSE]
}

cat("After goodSamplesGenes, datExpr0 dimensions:\n")
print(dim(datExpr0))

datExpr <- datExpr0

# -------------------------------------------------------------------
# 6) Sample clustering (outlier check) + save figure
# -------------------------------------------------------------------

sampleTree <- hclust(dist(datExpr), method = "average")

png("figures/GSE48301_WGCNA_sampleClustering.png",
    width = 2000, height = 1500, res = 300)
plot(sampleTree,
     main = "Sample clustering (GSE48301)",
     xlab = "", sub = "", cex = 0.7)
dev.off()

# For now keep all samples; if a clear outlier appears,
# I can drop it manually and re-run.

# -------------------------------------------------------------------
# 7) Choose soft-thresholding power + save figure
# -------------------------------------------------------------------

powers <- c(1:10, seq(12, 20, 2))
sft    <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

softPower <- sft$powerEstimate
if (is.na(softPower)) {
  softPower <- 6
  cat("pickSoftThreshold did not converge; falling back to softPower = 6.\n")
} else {
  cat("Chosen soft-thresholding power:", softPower, "\n")
}

png("figures/GSE48301_WGCNA_pickSoftThreshold.png",
    width = 2000, height = 1500, res = 300)
par(mfrow = c(1, 2))

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.8)
abline(h = 0.8, col = "red")

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, cex = 0.8)

par(mfrow = c(1, 1))
dev.off()

# -------------------------------------------------------------------
# 8) Network construction and dynamic tree cut
# -------------------------------------------------------------------

adjacency <- WGCNA::adjacency(datExpr, power = softPower)
TOM       <- WGCNA::TOMsimilarity(adjacency)
dissTOM   <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

png("figures/GSE48301_WGCNA_geneDendro_modules.png",
    width = 2000, height = 1500, res = 300)
WGCNA::plotDendroAndColors(
  geneTree,
  WGCNA::labels2colors(dynamicTreeCut::cutreeDynamic(
    dendro          = geneTree,
    distM           = dissTOM,
    deepSplit       = 3,      # slightly aggressive, to see more modules
    pamRespectsDend = FALSE,
    minClusterSize  = 15
  )),
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

# But I also need the colors object explicitly:

minModuleSize <- 15
dynamicMods <- dynamicTreeCut::cutreeDynamic(
  dendro          = geneTree,
  distM           = dissTOM,
  deepSplit       = 3,
  pamRespectsDend = FALSE,
  minClusterSize  = minModuleSize
)
dynamicColors <- WGCNA::labels2colors(dynamicMods)

cat("Module colours (including grey):\n")
print(table(dynamicColors))

# -------------------------------------------------------------------
# 9) Module eigengenes and merging
# -------------------------------------------------------------------

MEList <- WGCNA::moduleEigengenes(datExpr, colors = dynamicColors)
MEs    <- WGCNA::orderMEs(MEList$eigengenes)

moduleLevels <- setdiff(unique(dynamicColors), "grey")
cat("Number of initial modules (excluding grey):",
    length(moduleLevels), "\n")

if (length(moduleLevels) < 2) {
  # Only one real module; nothing meaningful to merge
  cat("Only one real module detected; skipping module merging.\n")
  moduleColors <- dynamicColors
  
} else {
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")
  
  png("figures/GSE48301_WGCNA_MEs_clustering.png",
      width = 2000, height = 1500, res = 300)
  plot(METree,
       main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h = 0.15, col = "red")
  dev.off()
  
  merge <- WGCNA::mergeCloseModules(
    datExpr,
    dynamicColors,
    cutHeight = 0.15,
    verbose  = 3
  )
  
  mergedColors <- merge$colors
  mergedMEs    <- WGCNA::orderMEs(merge$newMEs)
  
  cat("Number of modules before merging:",
      length(table(dynamicColors)), "\n")
  cat("Number of modules after  merging:",
      length(table(mergedColors)), "\n")
  
  moduleColors <- mergedColors
  MEs          <- mergedMEs
}

# -------------------------------------------------------------------
# 10) Module–trait relationships + heatmap figure
# -------------------------------------------------------------------

traitData <- traitData[rownames(datExpr), , drop = FALSE]

moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitP   <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(datExpr))

cat("Module–trait correlation matrix:\n")
print(moduleTraitCor)

textMatrix <- paste0(
  signif(moduleTraitCor, 2), "\n(",
  signif(moduleTraitP, 1), ")"
)
dim(textMatrix) <- dim(moduleTraitCor)

png("figures/GSE48301_WGCNA_moduleTrait_heatmap.png",
    width = 2000, height = 1500, res = 300)

WGCNA::labeledHeatmap(
  Matrix      = moduleTraitCor,
  xLabels     = colnames(traitData),
  yLabels     = rownames(moduleTraitCor),
  ySymbols    = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors      = WGCNA::blueWhiteRed(50),
  textMatrix  = textMatrix,
  main        = "Module–trait relationships (PCOS vs Control)"
)

dev.off()

# -------------------------------------------------------------------
# 11) Export results
# -------------------------------------------------------------------

geneModule_df <- data.frame(
  GeneID = colnames(datExpr),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

geneModule_file <- file.path("results", "GSE48301_WGCNA_gene_module_assignment.csv")
write.csv(geneModule_df, geneModule_file, row.names = FALSE)

moduleTrait_df <- data.frame(
  Module = rownames(moduleTraitCor),
  moduleTraitCor,
  moduleTraitP
)

moduleTrait_file <- file.path("results", "GSE48301_WGCNA_module_trait_correlations.csv")
write.csv(moduleTrait_df, moduleTrait_file, row.names = FALSE)

saveRDS(datExpr,      file.path("data", "GSE48301_WGCNA_datExpr.rds"))
saveRDS(moduleColors, file.path("data", "GSE48301_WGCNA_moduleColors.rds"))
saveRDS(MEs,          file.path("data", "GSE48301_WGCNA_MEs.rds"))

message("Saved gene–module assignment to: ", geneModule_file)
message("Saved module–trait correlations to: ", moduleTrait_file)
message("Saved WGCNA expression and module objects in 'data/' directory.")
message("Figures written to 'figures/' directory.")

# End of 04_WGCNA_network_and_modules_GSE48301.R
