# =============================================================================
# 05_hdwgcna_gsea.R
#
# Purpose: Gene Set Enrichment Analysis on hdWGCNA differentially expressed
#          modules (DMEs) using mouse-native MSigDB GO gene sets
# Input:   output/hdwgcna_results.rds
# Output:  output/gsea_results.rds
#          figures/gsea_barplot_Oligodendrocytes.pdf
#          figures/gsea_barplot_Glutamatergic_Neurons.pdf
#          figures/gsea_barplot_Gabaergic_Neurons.pdf
#          figures/gsea_barplot_Oligodendrocyte_Progenitor_Cells.pdf
#          figures/gsea_barplot_Endothelial_Cells.pdf
#          figures/gsea_barplot_Astrocytes.pdf
#
# Module selection criteria:
#   - Adjusted p-value < 0.05 (DME threshold)
#   - All significant modules included for completeness
#   - Top 8 GO terms per category ranked by absolute NES shown in figures
#
# Author: Karthikeya Kodali
# =============================================================================

library(Seurat)
library(hdWGCNA)
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)

# =============================================================================
# 1. LOAD DATA AND DME RESULTS
# =============================================================================

data <- LoadSeuratRds("output/hdwgcna_results.rds")

# load DME results - assumed to be computed in hdWGCNA pipeline
# DMEs_all should contain columns: celltype, module, avg_log2FC, p_val_adj
# if not already in environment, recompute:
cell_types_keep <- c(
  "Oligodendrocytes",
  "Glutamatergic Neurons",
  "Astrocytes",
  "Endothelial Cells",
  "Gabaergic Neurons",
  "Microglia"
)

DMEs_all <- data.frame()
for (ct in cell_types_keep) {
  group1 <- data@meta.data %>%
    subset(celltype == ct & condition == "WT_SD") %>% rownames
  group2 <- data@meta.data %>%
    subset(celltype == ct & condition == "WT_SD_Ctrl") %>% rownames
  
  if (length(group1) == 0 | length(group2) == 0) next
  
  cur_DMEs <- FindDMEs(
    data,
    barcodes1    = group1,
    barcodes2    = group2,
    test.use     = "wilcox",
    pseudocount.use = 0.01,
    wgcna_name   = ct
  )
  cur_DMEs$celltype <- ct
  DMEs_all <- rbind(DMEs_all, cur_DMEs)
}

# =============================================================================
# 2. LOAD MOUSE-NATIVE GO GENE SETS
# =============================================================================

# using mouse-native MSigDB (MM) to avoid ortholog mapping
go_gene_sets_mouse <- msigdbr(
  db_species = "MM",
  species    = "Mus musculus",
  collection = "M5"
) %>%
  filter(gs_subcollection %in% c("GO:BP", "GO:CC", "GO:MF")) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# verify overlap with dataset genes
all_mods    <- GetModules(data, wgcna_name = "Oligodendrocytes")
module_genes <- all_mods$gene_name
message("Gene overlap with MSigDB: ", 
        sum(module_genes %in% unique(unlist(go_gene_sets_mouse))),
        " / ", length(module_genes))
# expected: ~8817 genes overlap

# =============================================================================
# 3. RUN GSEA ON ALL SIGNIFICANT DME MODULES
# =============================================================================

# get all significant modules
sig_mods <- DMEs_all %>%
  filter(p_val_adj < 0.05) %>%
  select(celltype, module)

message("Running GSEA on ", nrow(sig_mods), " significant modules")

gsea_results <- list()

for (i in 1:nrow(sig_mods)) {
  ct  <- sig_mods$celltype[i]
  mod <- sig_mods$module[i]
  
  # get all genes in network for this cell type
  all_mods <- GetModules(data, wgcna_name = ct)
  
  # rank genes by kME for this module
  kme_col <- paste0("kME_", mod)
  
  ranked_genes <- all_mods %>%
    select(gene_name, all_of(kme_col)) %>%
    arrange(desc(.data[[kme_col]])) %>%
    deframe()
  
  # run fgsea with multilevel p-value estimation
  res <- fgseaMultilevel(
    pathways = go_gene_sets_mouse,
    stats    = ranked_genes,
    minSize  = 10,
    maxSize  = 500
  )
  
  gsea_results[[paste(ct, mod, sep = "_")]] <- res %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
  message("Done: ", ct, " - ", mod, " (",
          nrow(gsea_results[[paste(ct, mod, sep = "_")]]),
          " significant pathways)")
}

# save results
saveRDS(gsea_results, "output/gsea_results.rds")

# print summary
for (name in names(gsea_results)) {
  cat(name, ":", nrow(gsea_results[[name]]), "significant pathways\n")
}

# =============================================================================
# 4. PREPARE PLOTTING FUNCTIONS
# =============================================================================

# all significant modules organized by cell type
all_modules_list <- list(
  "Oligodendrocytes" = c(
    "Oligodendrocytes_Oligodendrocytes-M7",
    "Oligodendrocytes_Oligodendrocytes-M24",
    "Oligodendrocytes_Oligodendrocytes-M25",
    "Oligodendrocytes_Oligodendrocytes-M3",
    "Oligodendrocytes_Oligodendrocytes-M4",
    "Oligodendrocytes_Oligodendrocytes-M16"
  ),
  "Glutamatergic Neurons" = c(
    "Glutamatergic Neurons_GlutamatergicNeurons-M1",
    "Glutamatergic Neurons_GlutamatergicNeurons-M2",
    "Glutamatergic Neurons_GlutamatergicNeurons-M3",
    "Glutamatergic Neurons_GlutamatergicNeurons-M4",
    "Glutamatergic Neurons_GlutamatergicNeurons-M5"
  ),
  "Gabaergic Neurons" = c(
    "Gabaergic Neurons_GabaergicNeurons-M2",
    "Gabaergic Neurons_GabaergicNeurons-M9",
    "Gabaergic Neurons_GabaergicNeurons-M11",
    "Gabaergic Neurons_GabaergicNeurons-M17",
    "Gabaergic Neurons_GabaergicNeurons-M18"
  ),
  "Oligodendrocyte Progenitor Cells" = c(
    "Oligodendrocyte Progenitor Cells_OligodendrocyteProgenitorCells-M5",
    "Oligodendrocyte Progenitor Cells_OligodendrocyteProgenitorCells-M26",
    "Oligodendrocyte Progenitor Cells_OligodendrocyteProgenitorCells-M30"
  ),
  "Endothelial Cells" = c(
    "Endothelial Cells_EndothelialCells-M1",
    "Endothelial Cells_EndothelialCells-M3"
  ),
  "Astrocytes" = c(
    "Astrocytes_Astrocytes-M1"
  )
)

# function to prepare GSEA data for plotting
# selects top n_terms by absolute NES per GO category per module
prepare_gsea_plot_data <- function(module_keys, n_terms = 8) {
  
  do.call(rbind, lapply(module_keys, function(key) {
    res <- gsea_results[[key]]
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    mod_name <- gsub(".*_", "", key)
    
    # assign GO category
    res$category <- case_when(
      grepl("^GOBP_", res$pathway) ~ "Biological Process",
      grepl("^GOCC_", res$pathway) ~ "Cellular Component",
      grepl("^GOMF_", res$pathway) ~ "Molecular Function",
      TRUE ~ "Other"
    )
    
    # clean pathway names
    res$pathway_clean <- res$pathway %>%
      gsub("^GOBP_|^GOCC_|^GOMF_", "", .) %>%
      gsub("_", " ", .) %>%
      str_to_title() %>%
      str_wrap(width = 25)
    
    # top n_terms by absolute NES per category
    res %>%
      filter(category != "Other") %>%
      group_by(category) %>%
      slice_max(order_by = abs(NES), n = n_terms) %>%
      ungroup() %>%
      mutate(module = mod_name)
  }))
}

# function to generate GSEA barplot for one cell type
make_gsea_barplot <- function(ct, module_keys, n_terms = 8) {
  
  plot_data <- prepare_gsea_plot_data(module_keys, n_terms)
  if (is.null(plot_data) || nrow(plot_data) == 0) return(NULL)
  
  plot_data <- plot_data %>%
    arrange(category, module, NES) %>%
    mutate(
      pathway_clean = factor(pathway_clean, levels = unique(pathway_clean)),
      direction     = ifelse(NES > 0, "Upregulated in SD", "Downregulated in SD")
    )
  
  ggplot(plot_data, aes(x = NES, y = pathway_clean, fill = direction)) +
    geom_col(aes(alpha = -log10(padj))) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    scale_fill_manual(
      values = c(
        "Upregulated in SD"   = "#C0392B",
        "Downregulated in SD" = "#2980B9"
      ),
      name = ""
    ) +
    scale_alpha_continuous(
      name  = "-log10(padj)",
      range = c(0.3, 1)
    ) +
    facet_grid(category ~ module, scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(
      title = ct,
      x     = "Normalized Enrichment Score (NES)",
      y     = ""
    ) +
    theme(
      axis.text.y      = element_text(size = 9,  color = "black", face = "bold"),
      axis.text.x      = element_text(size = 10, color = "black", face = "bold"),
      axis.title.x     = element_text(size = 13, face = "bold"),
      strip.text.x     = element_text(size = 12, face = "bold"),
      strip.text.y     = element_text(size = 11, face = "bold", angle = 0),
      plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.title     = element_text(size = 12, face = "bold"),
      legend.text      = element_text(size = 11, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(0.5, "cm")
    )
}

# =============================================================================
# 5. GENERATE GSEA BARPLOTS PER CELL TYPE
# =============================================================================
# Top 8 GO terms per category ranked by absolute NES
# NES > 0 = upregulated in SD, NES < 0 = downregulated in SD

for (ct in names(all_modules_list)) {
  p <- make_gsea_barplot(ct, all_modules_list[[ct]], n_terms = 8)
  if (!is.null(p)) {
    fname  <- paste0("figures/gsea_barplot_", gsub(" ", "_", ct), ".pdf")
    n_mods <- length(all_modules_list[[ct]])
    width  <- 5 + (n_mods * 4)
    height <- 4 + (8 * 3 * 0.6)
    
    pdf(fname, width = width, height = height)
    print(p)
    dev.off()
    message("Saved: ", fname)
  }
}

message("GSEA analysis complete.")