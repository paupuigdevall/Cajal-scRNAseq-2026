# Environment setup for Practical 6: Cell Type Prediction
# Single-Cell and Spatial Omics Course 2026
#
# Verified against Noppe session (R 4.4.2, Bioconductor 3.20, Ubuntu 24.04)
# sessionInfo collected: 2026-03

# --- Bioconductor (version 3.20 matches R 4.4.x) ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install(c(
  "SingleR",              # 2.8.0
  "celldex",              # 1.16.0
  "BiocParallel",         # 1.40.2
  "SummarizedExperiment", # 1.36.0
  "SingleCellExperiment", # 1.28.1
  "GenomicRanges",        # 1.58.0
  "IRanges",              # 2.40.1
  "S4Vectors",            # 0.44.0
  "BiocGenerics",         # 0.52.0
  "Biobase",              # 2.66.0
  "pwalign",              # required by Azimuth
  "TFBSTools"             # required by Azimuth
), update = FALSE)

# --- CRAN ---
install.packages(c(
  "Seurat",       # 5.2.1
  "SeuratObject", # 5.0.2
  "dplyr",        # 1.1.4
  "ggplot2",      # 3.5.1
  "patchwork",    # 1.3.0
  "pheatmap",     # 1.0.12
  "fgsea",        # 1.32.4  (also on Bioconductor)
  "HGNChelper",   # 0.8.15
  "remotes",      # 2.5.0
  "sp",           # 2.2-0
  "matrixStats"   # 1.5.0
))

# --- GitHub ---
remotes::install_github("immunogenomics/presto@v1.0.0")          # 1.0.0
remotes::install_github("powellgenomicslab/scPred")               # 1.9.2
remotes::install_github("satijalab/azimuth")                      # latest
remotes::install_github("satijalab/seurat-data")                  # latest
