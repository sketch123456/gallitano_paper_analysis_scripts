# 01_preprocessing.R
# 
# Purpose: Convert processed Seurat object to loom format for pySCENIC
# Input: Seurat object (seu) - provided by Gallitano lab
# Output: .loom file for use in 02_pySCENIC.ipynb
#         "output/NS_only.rds"
#         "output/SD_only.rds"
#         .loom file of only SD cells for use in 02_pySCENIC.ipynb
#         .loom file of only NS cells for use in 02_pySCENIC.ipynb
# Note: Raw data preprocessing (QC, normalization) was performed 
# upstream by the Gallitano lab prior to this analysis
#
# Author: Karthikeya Kodali

library(Seurat)
library(SeuratExtend)
library(reticulate)

# Load seurat object
data <- readRDS("data/seurat_object.rds")

#replace with cell types as needed
new.cluster.ids <- c(
  "0"  = "Endothelial Cells", "10" = "Endothelial Cells",
  "18" = "Endothelial Cells", "19" = "Endothelial Cells",
  "1"  = "Microglia", "11" = "Microglia", "22" = "Microglia",
  "2"  = "Astrocytes", "4"  = "Astrocytes", "7"  = "Astrocytes",
  "8"  = "Astrocytes", "23" = "Astrocytes",
  "3"  = "Glutamatergic Neurons", "5"  = "Glutamatergic Neurons",
  "6"  = "Glutamatergic Neurons", "13" = "Glutamatergic Neurons",
  "14" = "Glutamatergic Neurons", "17" = "Glutamatergic Neurons",
  "9"  = "Oligodendrocyte Progenitor Cells",
  "12" = "Oligodendrocytes", "21" = "Oligodendrocytes",
  "15" = "Gabaergic Neurons", "16" = "Gabaergic Neurons",
  "20" = "Gabaergic Neurons"
)

#setting labels of seurat object
data <- RenameIdents(data, new.cluster.ids)
data$celltype <- Idents(data)
data$condition_label <- ifelse(data$condition == "WT_SD", "Sleep Deprived", "Control")
data$celltype_condition <- paste(data$celltype, data$condition, sep = "_")

#Split conditions into separate Seurat objects
NS <- subset(data, data@meta.data$condition == "WT_SD_Ctrl")
SD <- subset(data, data@meta.data$condition != "WT_SD_Ctrl")

#sanity check
isTrue(length(colnames(NS)) + length(colnames(SD)) == length(colnames(data)))

# Convert to loom format for pySCENIC input
Seu2Loom(
  data,
  filename = "data/scenic_input.loom",
  add.normdata = FALSE,
  add.metadata = TRUE,
  layers = NULL,
  overwrite = FALSE
)

Seu2Loom(
  SD,
  filename = "data/SD_scenic_input.loom",
  add.normdata = FALSE,
  add.metadata = TRUE,
  layers = NULL,
  overwrite = FALSE
)

Seu2Loom(
  NS,
  filename = "data/NS_scenic_input.loom",
  add.normdata = FALSE,
  add.metadata = TRUE,
  layers = NULL,
  overwrite = FALSE
)

SaveSeuratRds(SD, "output/SD_only.rds")
SaveSeuratRds(NS, "output/NS_only.rds")
