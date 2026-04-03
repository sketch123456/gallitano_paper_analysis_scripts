# =============================================================================
# 04_cellchat.R
#
# Purpose: Cell-cell communication analysis using CellChat
# Input:   output/hdwgcna_results.rds
# Output:  output/cellchat_objectlist.rds
#          output/cellchat_combined.rds
#          output/cellchat_net_up.rds
#          output/cellchat_net_down.rds
#          figures/cellchat_compare_interactions.pdf
#          figures/cellchat_diff_circle.pdf
#          figures/cellchat_diff_heatmap.pdf
#          figures/cellchat_ranknet.pdf
#          figures/cellchat_signalingrole_scatter.pdf
#          figures/cellchat_outgoing_heatmap.pdf
#          figures/cellchat_incoming_heatmap.pdf
#          figures/cellchat_PDGF_circle.pdf
#          figures/cellchat_bubble_DEG.pdf
#          figures/cellchat_functional_embedding.pdf
#          figures/cellchat_pathway_similarity.pdf
#
# Author: Karthikeya Kodali
# =============================================================================

library(CellChat)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)

# =============================================================================
# 1. LOAD DATA AND CREATE CONDITION-SPECIFIC SEURAT OBJECTS
# =============================================================================

data <- LoadSeuratRds("output/hdwgcna_results.rds")

# subset by condition
SD_seu   <- subset(data, condition == "WT_SD")
Ctrl_seu <- subset(data, condition == "WT_SD_Ctrl")

# set cell identities
Idents(SD_seu)   <- "celltype"
Idents(Ctrl_seu) <- "celltype"

# add samples column from biological sample ID
SD_seu$samples   <- SD_seu$org.sample.id
Ctrl_seu$samples <- Ctrl_seu$org.sample.id

# =============================================================================
# 2. CREATE CELLCHAT OBJECTS
# =============================================================================

cellchat_SD <- createCellChat(
  object   = SD_seu,
  group.by = "celltype",
  assay    = "RNA"
)

cellchat_Ctrl <- createCellChat(
  object   = Ctrl_seu,
  group.by = "celltype",
  assay    = "RNA"
)

# =============================================================================
# 3. SET DATABASE - USE ALL CATEGORIES
# =============================================================================

CellChatDB.use <- CellChatDB.mouse

cellchat_SD@DB   <- CellChatDB.use
cellchat_Ctrl@DB <- CellChatDB.use

# =============================================================================
# 4. RUN INFERENCE PIPELINE ON EACH CONDITION
# =============================================================================

run_cellchat_pipeline <- function(cellchat) {
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

cellchat_SD   <- run_cellchat_pipeline(cellchat_SD)
cellchat_Ctrl <- run_cellchat_pipeline(cellchat_Ctrl)

# =============================================================================
# 5. MERGE FOR COMPARISON
# =============================================================================

# SD first for internal indexing (SD = dataset 1, Control = dataset 2)
object.list <- list(
  SD      = cellchat_SD,
  Control = cellchat_Ctrl
)

cellchat_combined <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

# compute centrality scores for signaling role analysis
object.list <- lapply(object.list, function(x) {
  netAnalysis_computeCentrality(x, slot.name = "netP")
})

cellchat_SD   <- object.list[["SD"]]
cellchat_Ctrl <- object.list[["Control"]]

# =============================================================================
# 6. DEG-MAPPED L-R PAIR IDENTIFICATION
# =============================================================================

pos.dataset   <- "SD"
features.name <- paste0(pos.dataset, ".merged")

cellchat_combined <- identifyOverExpressedGenes(
  cellchat_combined,
  group.dataset     = "datasets",
  pos.dataset       = pos.dataset,
  features.name     = features.name,
  only.pos          = FALSE,
  thresh.pc         = 0.1,
  thresh.fc         = 0.05,
  thresh.p          = 0.05,
  group.DE.combined = FALSE
)

net      <- netMappingDEG(cellchat_combined, features.name = features.name,
                          variable.all = TRUE)
net.up   <- subsetCommunication(cellchat_combined, net = net,
                                datasets = "SD",
                                ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat_combined, net = net,
                                datasets = "Control",
                                ligand.logFC = -0.05, receptor.logFC = NULL)

# =============================================================================
# 7. SAVE OUTPUTS
# =============================================================================

saveRDS(object.list,       "output/cellchat_objectlist.rds")
saveRDS(cellchat_combined, "output/cellchat_combined.rds")
saveRDS(net.up,            "output/cellchat_net_up.rds")
saveRDS(net.down,          "output/cellchat_net_down.rds")

# =============================================================================
# 8. VISUALIZATIONS
# =============================================================================

# create swapped object list for differential plots
# (Control first so red = increased in SD)
object.list_plot <- list(
  Control = cellchat_Ctrl,
  SD      = cellchat_SD
)
cellchat_plot <- mergeCellChat(
  object.list_plot,
  add.names = names(object.list_plot)
)

# ── 8a. Compare total interactions ───────────────────────────────────────────
pdf("figures/cellchat_compare_interactions.pdf", width = 8, height = 5)
gg1 <- compareInteractions(
  cellchat_combined,
  show.legend = FALSE,
  group       = c(1, 2),
  color.use   = c("SD" = "firebrick", "Control" = "steelblue")
)
gg2 <- compareInteractions(
  cellchat_combined,
  show.legend = FALSE,
  group       = c(1, 2),
  measure     = "weight",
  color.use   = c("SD" = "firebrick", "Control" = "steelblue")
)
print(gg1 + gg2)
dev.off()

# ── 8b. Differential circle plots (red = increased in SD) ────────────────────
pdf("figures/cellchat_diff_circle.pdf", width = 14, height = 6)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(
  cellchat_plot,
  weight.scale = TRUE,
  label.edge   = FALSE,
  title.name   = "Differential number of interactions"
)
netVisual_diffInteraction(
  cellchat_plot,
  weight.scale = TRUE,
  measure      = "weight",
  label.edge   = FALSE,
  title.name   = "Differential interaction strength"
)
dev.off()

# ── 8c. Differential heatmap (red = increased in SD) ─────────────────────────
gg3 <- netVisual_heatmap(
  cellchat_plot,
  title.name = "Differential number of interactions"
)
gg4 <- netVisual_heatmap(
  cellchat_plot,
  measure    = "weight",
  title.name = "Differential interaction strength"
)
pdf("figures/cellchat_diff_heatmap.pdf", width = 14, height = 6)
ComplexHeatmap::draw(gg3 + gg4, ht_gap = unit(0.5, "cm"))
dev.off()

# ── 8d. rankNet comparison ────────────────────────────────────────────────────
pdf("figures/cellchat_ranknet.pdf", width = 14, height = 8)
gg5 <- rankNet(
  cellchat_combined,
  mode      = "comparison",
  measure   = "weight",
  stacked   = TRUE,
  do.stat   = TRUE,
  color.use = c("SD" = "firebrick", "Control" = "steelblue")
)
gg6 <- rankNet(
  cellchat_combined,
  mode      = "comparison",
  measure   = "weight",
  stacked   = FALSE,
  do.stat   = TRUE,
  color.use = c("SD" = "firebrick", "Control" = "steelblue")
)
print(gg5 + gg6)
dev.off()

# ── 8e. Signaling role scatter ────────────────────────────────────────────────
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))

gg_scatter <- list()
for (i in 1:length(object.list)) {
  gg_scatter[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title         = names(object.list)[i],
    weight.MinMax = weight.MinMax
  )
}

pdf("figures/cellchat_signalingrole_scatter.pdf", width = 14, height = 6)
patchwork::wrap_plots(plots = gg_scatter)
dev.off()

# ── 8f. Outgoing and incoming signaling role heatmaps ────────────────────────
pathway.union <- union(
  object.list[[1]]@netP$pathways,
  object.list[[2]]@netP$pathways
)

# outgoing
ht_out1 <- netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "outgoing",
  signaling = pathway.union,
  title     = names(object.list)[1],
  width = 6, height = 12
)
ht_out2 <- netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "outgoing",
  signaling = pathway.union,
  title     = names(object.list)[2],
  width = 6, height = 12
)
pdf("figures/cellchat_outgoing_heatmap.pdf", width = 14, height = 16)
ComplexHeatmap::draw(ht_out1 + ht_out2, ht_gap = unit(0.5, "cm"))
dev.off()

# incoming
ht_in1 <- netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "incoming",
  signaling     = pathway.union,
  title         = names(object.list)[1],
  width = 6, height = 12,
  color.heatmap = "GnBu"
)
ht_in2 <- netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "incoming",
  signaling     = pathway.union,
  title         = names(object.list)[2],
  width = 6, height = 12,
  color.heatmap = "GnBu"
)
pdf("figures/cellchat_incoming_heatmap.pdf", width = 14, height = 16)
ComplexHeatmap::draw(ht_in1 + ht_in2, ht_gap = unit(0.5, "cm"))
dev.off()

# ── 8g. PDGF pathway circle plot ─────────────────────────────────────────────
weight.max.pdgf <- getMaxWeight(
  object.list,
  slot.name = c("netP"),
  attribute = "PDGF"
)

pdf("figures/cellchat_PDGF_circle.pdf", width = 12, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(
    object.list[[i]],
    signaling       = "PDGF",
    layout          = "circle",
    edge.weight.max = weight.max.pdgf[1],
    edge.width.max  = 10,
    signaling.name  = paste("PDGF", names(object.list)[i])
  )
}
dev.off()

# ── 8h. DEG-mapped bubble plots ───────────────────────────────────────────────
pairLR.use.up   <- net.up[,   "interaction_name", drop = FALSE]
pairLR.use.down <- net.down[, "interaction_name", drop = FALSE]

gg_bub1 <- netVisual_bubble(
  cellchat_combined,
  pairLR.use     = pairLR.use.up,
  comparison     = c(1, 2),
  angle.x        = 90,
  remove.isolate = TRUE,
  title.name     = "Upregulated signaling in SD"
)

gg_bub2 <- netVisual_bubble(
  cellchat_combined,
  pairLR.use     = pairLR.use.down,
  comparison     = c(1, 2),
  angle.x        = 90,
  remove.isolate = TRUE,
  title.name     = "Downregulated signaling in SD"
)

pdf("figures/cellchat_bubble_DEG.pdf", width = 20, height = 14)
print(gg_bub1 + gg_bub2)
dev.off()

# ── 8i. Functional similarity embedding ───────────────────────────────────────
cellchat_combined <- computeNetSimilarityPairwise(
  cellchat_combined, type = "functional"
)
cellchat_combined <- netEmbedding(
  cellchat_combined, type = "functional"
)
cellchat_combined <- netClustering(
  cellchat_combined, type = "functional"
)

pdf("figures/cellchat_functional_embedding.pdf", width = 8, height = 6)
netVisual_embeddingPairwise(
  cellchat_combined, type = "functional", label.size = 3.5
)
dev.off()

pdf("figures/cellchat_pathway_similarity.pdf", width = 8, height = 5)
rankSimilarity(cellchat_combined, type = "functional")
dev.off()

# resave combined object with embedding
saveRDS(cellchat_combined, "output/cellchat_combined.rds")