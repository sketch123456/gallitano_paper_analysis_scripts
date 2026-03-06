# 02_hdWGCNA.R
#
# Purpose: Run hdWGCNA co-expression network analysis across all cell types
# Input: data/seurat_object.rds
# Output: output/hdwgcna_results.rds
#
# Performs metacell construction, soft power estimation, network construction,
# and module eigengene/connectivity calculation for each cell type iteratively.
# Failed cell types are skipped with error logging.
#
# Author: Karthikeya Kodali

library(Seurat)
library(SeuratExtend)
library(hdWGCNA)
library(WGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)

# Load processed Seurat object
data <- LoadSeuratRds("data/seurat_object.rds")

# Set default assay
DefaultAssay(data) <- "RNA"

# Get list of unique cell types
cell_types <- unique(data$celltype)

# Run hdWGCNA for each cell type
for (ct in cell_types) {
  message("\n========== Processing: ", ct, " ==========")
  
  tryCatch({
    # Setup hdWGCNA
    data <- SetupForWGCNA(
      data,
      gene_select = "fraction",
      fraction = 0.05,
      wgcna_name = "hdWGCNA"
    )
    
    # Create metacells per sample x celltype combination
    data <- MetacellsByGroups(
      seurat_obj = data,
      group.by = c("celltype", "org.sample.id"),
      reduction = "umap",
      k = 15,
      max_shared = 10,
      ident.group = "celltype",
      wgcna_name = ct
    )
    
    # Normalize metacells
    data <- NormalizeMetacells(data, wgcna_name = ct)
    
    # Set expression data
    data <- SetDatExpr(
      data,
      group_name = ct,
      group.by = "celltype",
      assay = "RNA",
      layer = "data",
      use_metacells = TRUE,
      wgcna_name = ct
    )
    
    # Determine optimal soft power threshold
    data <- TestSoftPowers(data, networkType = "signed", wgcna_name = ct)
    
    # Construct co-expression network
    data <- ConstructNetwork(
      data,
      tom_name = gsub(" ", "_", ct),
      overwrite_tom = TRUE,
      wgcna_name = ct
    )
    
    # Scale data for eigengene calculation
    data <- ScaleData(data, features = VariableFeatures(data))
    
    # Calculate module eigengenes and connectivity
    data <- ModuleEigengenes(
      data,
      group.by.vars = "org.sample.id",
      wgcna_name = ct
    )
    
    data <- ModuleConnectivity(
      data,
      group.by = "celltype",
      group_name = ct,
      wgcna_name = ct
    )
    
    # Rename modules with cell type prefix (e.g. Neuron-M1)
    prefix <- gsub(" ", "", ct)
    data <- ResetModuleNames(
      data,
      new_name = paste0(prefix, "-M"),
      wgcna_name = ct
    )
    
    message("Completed: ", ct)
    
  }, error = function(e) {
    message("ERROR in ", ct, ": ", e$message)
    message("Skipping and continuing...")
  })
}

# Save results
SaveSeuratRds(data, "output/hdwgcna_results.rds")