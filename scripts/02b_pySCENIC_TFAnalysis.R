# 05_pySCENIC_analysis.R
#
# Purpose: Import and analyze pySCENIC results within Seurat framework
# Input:   data/seurat_object.rds
#          data/pyscenic_integrated-output.loom
# Output:  figures/tf_zscore_heatmap.pdf
#          figures/TF_waterfall_[celltype].pdf
#          figures/jun_regulon_violin.pdf
#          output/regulons_and_targets_wide.csv
#
# Analyzes transcription factor regulon activity (AUCell scores) from pySCENIC
# comparing WT_SD vs WT_SD_Ctrl conditions across cell types.
# Cluster annotations, differential TF activity, and Jun regulon
# visualizations are included.
#
# Author: Karth Kodali
# Gallitano Lab - University of Arizona
# Date: 2025

library(Seurat)
library(SeuratExtend)
library(ggplot2)
library(ggrepel)
library(tidyverse)

# ── 1. Load data and import pySCENIC results ──────────────────────────────────

data <- LoadSeuratRds("data/seurat_object.rds")
data <- ImportPyscenicLoom("data/pyscenic_integrated-output.loom", seu = data)

# Extract AUC matrix (cells x TFs) and regulon gene lists
tf_auc       <- data@misc$SCENIC$RegulonsAUC
tf_gene_list <- data@misc$SCENIC$Regulons

# ── 2. Annotate clusters with cell type labels ────────────────────────────────
# Cluster assignments specific to this dataset (sleep deprivation study)

new.cluster.ids <- c(
  "0"  = "Endothelial Cells",
  "10" = "Endothelial Cells",
  "18" = "Endothelial Cells",
  "19" = "Endothelial Cells",
  "1"  = "Microglia",
  "11" = "Microglia",
  "22" = "Microglia",
  "2"  = "Astrocytes",
  "4"  = "Astrocytes",
  "7"  = "Astrocytes",
  "8"  = "Astrocytes",
  "23" = "Astrocytes",
  "3"  = "Glutamatergic Neurons",
  "5"  = "Glutamatergic Neurons",
  "6"  = "Glutamatergic Neurons",
  "13" = "Glutamatergic Neurons",
  "14" = "Glutamatergic Neurons",
  "17" = "Glutamatergic Neurons",
  "9"  = "Oligodendrocyte Progenitor Cells",
  "12" = "Oligodendrocytes",
  "21" = "Oligodendrocytes",
  "15" = "GABAergic Neurons",
  "16" = "GABAergic Neurons",
  "20" = "GABAergic Neurons"
)

data <- RenameIdents(data, new.cluster.ids)
data$celltype <- Idents(data)
Idents(data) <- "celltype"

# ── 3. TF z-score heatmap across clusters ────────────────────────────────────

tf_zscore <- CalcStats(
  tf_auc,
  f           = data$celltype,
  order       = "p",
  p.threshold = 0.1,
  n           = 10,
  t           = TRUE
)

pdf("figures/tf_zscore_heatmap.pdf", width = 12, height = 8)
Heatmap(tf_zscore, lab_fill = "zscore")
dev.off()

# ── 4. Differential TF activity: WT_SD vs WT_SD_Ctrl per cell type ───────────
# Bar plot of top 10 up in SD and top 10 up in Ctrl per cell type

# Build combined dataframe of mean AUC difference per TF per cell type
all_volcano <- do.call(rbind, lapply(unique(data$celltype), function(ct) {
  
  sd_cells   <- rownames(data@meta.data[data$celltype == ct & data$condition == "WT_SD", ])
  ctrl_cells <- rownames(data@meta.data[data$celltype == ct & data$condition == "WT_SD_Ctrl", ])
  
  auc_sd   <- tf_auc[sd_cells,   , drop = FALSE]
  auc_ctrl <- tf_auc[ctrl_cells, , drop = FALSE]
  
  avg_diff <- colMeans(auc_sd, na.rm = TRUE) - colMeans(auc_ctrl, na.rm = TRUE)
  
  data.frame(
    celltype = ct,
    TF       = names(avg_diff),
    avg_diff = as.numeric(avg_diff),
    stringsAsFactors = FALSE
  )
}))

# Plot top 20 differentially active TFs per cell type
for (ct in unique(all_volcano$celltype)) {
  df <- subset(all_volcano, celltype == ct)
  
  # Top 10 up in SD, top 10 up in Ctrl
  top_up   <- df[order(df$avg_diff, decreasing = TRUE),  ][1:10, ]
  top_down <- df[order(df$avg_diff, decreasing = FALSE), ][1:10, ]
  plot_df  <- rbind(top_up, top_down)
  
  plot_df$TF_clean  <- gsub("\\(\\+\\)", "", plot_df$TF)
  plot_df$direction <- ifelse(plot_df$avg_diff > 0, "Up in SD", "Up in Ctrl")
  plot_df           <- plot_df[order(plot_df$avg_diff), ]
  plot_df$TF_clean  <- factor(plot_df$TF_clean, levels = plot_df$TF_clean)
  
  p <- ggplot(plot_df, aes(avg_diff, TF_clean, fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Up in SD" = "firebrick", "Up in Ctrl" = "steelblue")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = ct,
      x     = "Mean AUC difference (SD - Ctrl)",
      y     = NULL,
      fill  = NULL
    ) +
    theme_classic()
  
  pdf(paste0("figures/TF_barplot_", gsub(" ", "_", ct), ".pdf"),
      width = 8, height = 6)
  print(p)
  dev.off()
}

# ── 5. Jun regulon activity violin plot ───────────────────────────────────────

celltype_order <- c(
  "Glutamatergic Neurons", "GABAergic Neurons",
  "Astrocytes",            "Microglia",
  "Oligodendrocytes",      "Oligodendrocyte Progenitor Cells",
  "Endothelial Cells"
)

# Recode condition labels for clean plot display
data$condition_label <- ifelse(
  data$condition == "WT_SD", "Sleep Deprived", "Control"
)

# Add Jun AUC to metadata
data$tf_plot <- tf_auc[colnames(data), "Jun(+)"]

# Set factor order
data$celltype <- factor(data$celltype, levels = celltype_order)

p_jun <- VlnPlot(
  data,
  features  = "tf_plot",
  group.by  = "celltype",
  split.by  = "condition_label",
  pt.size   = 0,
  cols      = c("Sleep Deprived" = "firebrick", "Control" = "steelblue")
) +
  labs(
    title = "Jun Regulon Activity",
    y     = "Regulon AUC",
    x     = NULL,
    fill  = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    plot.title   = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )

pdf("figures/jun_regulon_violin.pdf", width = 10, height = 6)
print(p_jun)
dev.off()

# ── 6. Export regulon target gene lists ───────────────────────────────────────

max_targets    <- max(sapply(tf_gene_list, length))

regulon_wide <- do.call(rbind, lapply(names(tf_gene_list), function(tf) {
  genes         <- tf_gene_list[[tf]]
  length(genes) <- max_targets
  c(TF = tf, genes)
}))

regulon_wide_df <- as.data.frame(regulon_wide, stringsAsFactors = FALSE)
colnames(regulon_wide_df) <- c("TF", paste0("target_", 1:max_targets))

write.csv(regulon_wide_df, "output/regulons_and_targets_wide.csv",
          row.names = FALSE)

message("pySCENIC analysis complete.")