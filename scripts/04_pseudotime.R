# 03_pseudotime.R
#
# Purpose: Pseudotime analysis of OPC to Oligodendrocyte differentiation
# Input:   output/hdwgcna_results.rds
# Output:  figures/pseudotime_palantir_dimplot.pdf
#          figures/pseudotime_density_palantir.pdf
#          figures/gene_trends_palantir.pdf
#          figures/pseudotime_density_slingshot.pdf
#          figures/gene_trends_slingshot.pdf
#          output/odc_lineage.rds
#
#Performs Pseudotime trajectory analysis using both Palantir and Slingshot
#on the initial full Seurat object. Results showed that OPC-ODC was the 
#only biologically meaningful trajectory, which was then further analyzed
#for effects of SD on progression and on gene expression over pseudotime.
#
# Author: Karthikeya Kodali

library(Seurat)
library(SeuratExtend)
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)

# ── 1. INITIAL FULL-DATASET PALANTIR ─────────────────────────────────────────
# Exploratory analysis on full dataset — not used for final analysis
# Justifies subsetting to OPC-ODC lineage only

data <- LoadSeuratRds("output/hdwgcna_results.rds")

# Set conda environment for Python-dependent Palantir (must change file path here)
options(reticulate.conda_binary = "C:/Users/kkoda/anaconda3/Scripts/conda.exe")
activate_python()

# Run diffusion map on full dataset
data <- Palantir.RunDM(data, reduction = "pca")

data <- RunUMAP(
  data,
  reduction      = "ms",
  dims           = 1:ncol(data@reductions$ms@cell.embeddings),
  reduction.name = "umap_ms"
)

# Start cell: highest Pdgfra OPC (most immature)
opc_cells  <- WhichCells(data, idents = "Oligodendrocyte Progenitor Cells")
start_cell <- names(which.max(data[["RNA"]]$data["Pdgfra", opc_cells]))

data <- Palantir.Pseudotime(data, start_cell = start_cell)

ps <- data@misc$Palantir$Pseudotime
data@meta.data[, colnames(ps)] <- ps

pdf("figures/pseudotime_full_dataset_exploratory.pdf", width = 10, height = 5)
DimPlot2(data,
         features  = c("Pseudotime", "Entropy"),
         reduction = "umap",
         cols      = list(Pseudotime = "D", Entropy = "D")
)
dev.off()
# Result: entropy = 0 throughout, confirming full dataset is not appropriate
# for trajectory analysis due to 7 discrete disconnected clusters

# ── 2. FULL DATASET SLINGSHOT ─────────────────────────────────────────────────
# Exploratory — confirms subsetting rationale

sce_full <- as.SingleCellExperiment(data)
reducedDim(sce_full, "MS") <- data@reductions$ms@cell.embeddings

sce_full <- slingshot(
  sce_full,
  clusterLabels = "celltype",
  reducedDim    = "MS",
  start.clus    = "Oligodendrocyte Progenitor Cells"
)
# Result: 4 lineages detected, only lineage 4 (OPC->ODC) is biologically valid
# All other lineages connect OPCs to neurons/endothelial - not meaningful
# Further justification for subsetting to OPC-ODC lineage

# ── 3. SUBSET OPC-OLIGODENDROCYTE LINEAGE ────────────────────────────────────

odc_lineage <- subset(data, celltype %in% c(
  "Oligodendrocyte Progenitor Cells",
  "Oligodendrocytes"
))

# Cell counts per condition
# OPC:  307 SD, 495 control
# ODC:  150 SD, 381 control
# Note: SD depletion more severe in mature ODCs — suggests impaired maturation
print(table(odc_lineage$celltype, odc_lineage$condition))

# ── 4. PALANTIR ON OPC-ODC LINEAGE ───────────────────────────────────────────

# Rerun PCA on subset — subsetting changes expression landscape
odc_lineage <- RunPCA(odc_lineage, assay = "RNA")
odc_lineage <- Palantir.RunDM(odc_lineage, reduction = "pca")

odc_lineage <- RunUMAP(
  odc_lineage,
  reduction      = "ms",
  dims           = 1:ncol(odc_lineage@reductions$ms@cell.embeddings),
  reduction.name = "umap_ms"
)

# Start cell: highest Pdgfra OPC
opc_cells  <- WhichCells(odc_lineage, idents = "Oligodendrocyte Progenitor Cells")
start_cell <- names(which.max(odc_lineage[["RNA"]]$data["Pdgfra", opc_cells]))

# Terminal cell: highest Mbp ODC (most mature)
odc_cells     <- WhichCells(odc_lineage, idents = "Oligodendrocytes")
terminal_cell <- names(which.max(odc_lineage[["RNA"]]$data["Mbp", odc_cells]))

# Verify start and terminal cell locations
odc_lineage$is_start    <- colnames(odc_lineage) == start_cell
odc_lineage$is_terminal <- colnames(odc_lineage) == terminal_cell

pdf("figures/start_terminal_cell_locations.pdf", width = 12, height = 5)
DimPlot(odc_lineage, reduction = "umap_ms", group.by = "is_start",
        cols = c("grey80", "red"), order = TRUE) + ggtitle("Start cell") |
  DimPlot(odc_lineage, reduction = "umap_ms", group.by = "is_terminal",
          cols = c("grey80", "red"), order = TRUE) + ggtitle("Terminal cell")
dev.off()

# Run pseudotime
odc_lineage <- Palantir.Pseudotime(
  odc_lineage,
  start_cell      = start_cell,
  terminal_states = c("fate_ODC" = terminal_cell)
)

ps <- odc_lineage@misc$Palantir$Pseudotime
odc_lineage@meta.data[, colnames(ps)] <- ps

pdf("figures/pseudotime_palantir_dimplot.pdf", width = 10, height = 5)
DimPlot2(odc_lineage,
         features  = c("Pseudotime", "Entropy"),
         reduction = "umap_ms",
         cols      = list(Pseudotime = "D", Entropy = "D")
)
dev.off()

# ── 4a. Palantir: Compare conditions ─────────────────────────────────────────

odc_meta <- odc_lineage@meta.data %>% filter(!is.na(Pseudotime))

p_density_palantir <- ggplot(
  odc_meta, aes(x = Pseudotime, fill = condition, color = condition)
) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values  = c("WT_SD" = "firebrick", "WT_SD_Ctrl" = "steelblue")) +
  scale_color_manual(values = c("WT_SD" = "firebrick", "WT_SD_Ctrl" = "steelblue")) +
  theme_bw() +
  labs(
    title = "OPC-Oligodendrocyte Pseudotime: SD vs Control (Palantir)",
    x     = "Pseudotime",
    y     = "Density"
  )

pdf("figures/pseudotime_density_palantir.pdf", width = 8, height = 5)
print(p_density_palantir)
dev.off()

# Wilcoxon test: SD vs Ctrl pseudotime distribution
wilcox_palantir <- wilcox.test(
  odc_lineage$Pseudotime[odc_lineage$condition == "WT_SD"],
  odc_lineage$Pseudotime[odc_lineage$condition == "WT_SD_Ctrl"]
)
print(wilcox_palantir)
# Result: p = 9.448e-09
# SD cells enriched at early pseudotime (OPC state)
# Control cells enriched at late pseudotime (mature ODC state)

# ── 4b. Palantir: Gene trend curves ──────────────────────────────────────────

genes_of_interest <- c(
  "Pdgfra", "Cspg4",      # OPC markers
  "Mbp", "Plp1", "Mobp",  # mature ODC markers
  "Jun", "Egr1",           # IEGs / TFs of interest
  "Sox10", "Olig2"         # lineage markers
)

# Split by condition
odc_SD   <- subset(odc_lineage, condition == "WT_SD")
odc_Ctrl <- subset(odc_lineage, condition == "WT_SD_Ctrl")

# Copy pseudotime to subsets
odc_SD@misc$Palantir$Pseudotime   <- ps[rownames(odc_SD@meta.data), ]
odc_Ctrl@misc$Palantir$Pseudotime <- ps[rownames(odc_Ctrl@meta.data), ]

p_SD_pal <- GeneTrendCurve.Palantir(
  odc_SD,
  features        = genes_of_interest,
  pseudotime.data = "Pseudotime",
  method          = "gam",
  se              = TRUE
) + ggtitle("Sleep Deprived")

p_Ctrl_pal <- GeneTrendCurve.Palantir(
  odc_Ctrl,
  features        = genes_of_interest,
  pseudotime.data = "Pseudotime",
  method          = "gam",
  se              = TRUE
) + ggtitle("Control")

pdf("figures/gene_trends_palantir.pdf", width = 10, height = 10)
print(p_SD_pal / p_Ctrl_pal)
dev.off()

# ── 5. SLINGSHOT ON OPC-ODC LINEAGE (validation) ─────────────────────────────

odc_lineage <- RunSlingshot(
  odc_lineage,
  group.by   = "celltype",
  start.clus = "Oligodendrocyte Progenitor Cells"
)

sling_full <- odc_lineage@misc$slingshot$PCA$SlingPseudotime
odc_lineage@meta.data[, colnames(sling_full)] <- as.data.frame(sling_full)

pdf("figures/pseudotime_slingshot_dimplot.pdf", width = 8, height = 5)
DimPlot2(odc_lineage, features = colnames(sling_full), cols = "C", theme = NoAxes())
dev.off()

# ── 5a. Slingshot: Compare conditions ────────────────────────────────────────

odc_meta <- odc_lineage@meta.data %>% filter(!is.na(slingPseudotime_1))
odc_meta$pt_scaled <- (odc_meta$slingPseudotime_1 - min(odc_meta$slingPseudotime_1)) /
  (max(odc_meta$slingPseudotime_1) - min(odc_meta$slingPseudotime_1))

p_density_sling <- ggplot(
  odc_meta, aes(x = pt_scaled, fill = condition, color = condition)
) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values  = c("WT_SD" = "firebrick", "WT_SD_Ctrl" = "steelblue")) +
  scale_color_manual(values = c("WT_SD" = "firebrick", "WT_SD_Ctrl" = "steelblue")) +
  theme_bw() +
  labs(
    title = "Slingshot Pseudotime: SD vs Control",
    x     = "Pseudotime (scaled 0-1)",
    y     = "Density"
  )

pdf("figures/pseudotime_density_slingshot.pdf", width = 8, height = 5)
print(p_density_sling)
dev.off()

wilcox_sling <- wilcox.test(
  odc_lineage$slingPseudotime_1[odc_lineage$condition == "WT_SD"],
  odc_lineage$slingPseudotime_1[odc_lineage$condition == "WT_SD_Ctrl"]
)
print(wilcox_sling)
# Result: p = 0.0006177
# Replicates Palantir finding

# ── 5b. Slingshot: Gene trend curves ─────────────────────────────────────────

odc_SD   <- subset(odc_lineage, condition == "WT_SD")
odc_Ctrl <- subset(odc_lineage, condition == "WT_SD_Ctrl")

odc_SD@misc$slingshot   <- odc_lineage@misc$slingshot
odc_Ctrl@misc$slingshot <- odc_lineage@misc$slingshot

odc_SD@misc$slingshot$PCA$SlingPseudotime   <- sling_full[rownames(odc_SD@meta.data),   , drop = FALSE]
odc_Ctrl@misc$slingshot$PCA$SlingPseudotime <- sling_full[rownames(odc_Ctrl@meta.data), , drop = FALSE]

p_SD_sling <- GeneTrendCurve.Slingshot(
  odc_SD, features = genes_of_interest
) + ggtitle("Sleep Deprived")

p_Ctrl_sling <- GeneTrendCurve.Slingshot(
  odc_Ctrl, features = genes_of_interest
) + ggtitle("Control")

pdf("figures/gene_trends_slingshot.pdf", width = 10, height = 10)
print(p_SD_sling / p_Ctrl_sling)
dev.off()

# Save lineage object for downstream use
SaveSeuratRds(odc_lineage, "output/odc_lineage.rds")