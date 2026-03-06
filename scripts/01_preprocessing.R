# 01_preprocessing.R
# 
# Purpose: Convert processed Seurat object to loom format for pySCENIC
# Input: Seurat object (seu) - provided by Gallitano lab
# Output: .loom file for use in 01_pySCENIC.ipynb
# Note: Raw data preprocessing (QC, normalization) was performed 
# upstream by the Gallitano lab prior to this analysis
#
# Author: Karthikeya Kodali

library(Seurat)
library(SeuratExtend)
library(reticulate)

# Load seurat object
seu <- readRDS("data/seurat_object.rds")

# Convert to loom format for pySCENIC input
Seu2Loom(
  seu,
  filename = "data/scenic_input.loom",
  add.normdata = FALSE,
  add.metadata = TRUE,
  layers = NULL,
  overwrite = FALSE
)