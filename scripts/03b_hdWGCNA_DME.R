# 02b_hdWGCNA_DME.R
#
# Purpose: Differential Module Eigengene (DME) analysis and visualization
# Input: output/hdwgcna_results.rds
# Output: figures/DME_volcano_[celltype].pdf
#         figures/hubgene_network_[celltype].pdf
#
# Identifies sex-differentiated co-expression modules per cell type,
# generates volcano plots of DMEs, and visualizes hub gene networks
# for significant modules.
#
# Author: Karthikeya Kodali

library(Seurat)
library(hdWGCNA)
library(tidyverse)
library(patchwork)

# Load hdWGCNA results
seurat_obj <- LoadSeuratRds("output/hdwgcna_results.rds")

# Get cell types
cell_types <- unique(seurat_obj$celltype)

for (ct in cell_types) {
  message("\n========== DME Analysis: ", ct, " ==========")
  
  tryCatch({
    
    # Define groups by condition
    # group1 = Sleep Deprived (WT_SD), group2 = Control (WT_SD_Ctrl)
    group1 <- seurat_obj@meta.data %>%
      subset(celltype == ct & condition == "WT_SD") %>%
      rownames()
    
    group2 <- seurat_obj@meta.data %>%
      subset(celltype == ct & condition == "WT_SD_Ctrl") %>%
      rownames()
    
    # Run DME analysis (Wilcoxon test)
    DMEs <- FindDMEs(
      seurat_obj,
      barcodes1 = group1,
      barcodes2 = group2,
      test.use = "wilcox",
      wgcna_name = ct
    )
    
    # Generate and save volcano plot
    p <- PlotDMEsVolcano(
      seurat_obj,
      DMEs,
      wgcna_name = ct
    )
    
    pdf(paste0("figures/DME_volcano_", gsub(" ", "_", ct), ".pdf"))
    print(p)
    dev.off()
    
    # Filter for significant DMEs
    # Criteria: adjusted p-value < 0.05 and absolute fold change > 0.25
    sig_DMEs <- DMEs %>%
      filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)
    
    if (nrow(sig_DMEs) == 0) {
      message("No significant DMEs for ", ct, ", skipping hub gene plot")
      next
    }
    
    # Get significant module names
    sig_mods <- unique(sig_DMEs$module)
    
    message("Significant modules for ", ct, ": ", paste(sig_mods, collapse = ", "))
    
    # Generate hub gene network plot for significant modules
    pdf(paste0("figures/hubgene_network_", gsub(" ", "_", ct), ".pdf"),
        width = 10, height = 10)
    
    HubGeneNetworkPlot(
      seurat_obj,
      n_hubs = 8,
      n_other = 2,
      edge_prop = 0.75,
      mods = sig_mods,
      wgcna_name = ct
    )
    
    dev.off()
    
    message("Completed DME analysis: ", ct)
    
  }, error = function(e) {
    message("ERROR in ", ct, ": ", e$message)
    message("Skipping and continuing...")
  })
}