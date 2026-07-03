###############################################################################
# Monocle3 拟时序 — 完整修正版
# 微胶质细胞子集, 合并分析, use_partition=TRUE, 手动root, 全局归一化
###############################################################################
library(Seurat)
library(monocle3)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
set.seed(42)

# ============================================================
# STEP 0-2: Extract & merge microglia, re-process from scratch
# ============================================================
cat("=== STEP 0-2: Build microglia Seurat object ===\n")
seu_cf   <- readRDS("data/processed/seu_cf_qc.rds")
seu_sham <- readRDS("data/processed/seu_sham_qc.rds")
seu_anno <- readRDS("data/processed/seu_mouse_annotated.rds")

strip_prefix <- function(x) sub("^(CF_|Sham_)", "", x)
mg_cf_barcodes   <- colnames(seu_anno)[seu_anno$cell_type == "Microglia" & seu_anno$orig.ident == "Mouse_CF"]
mg_sham_barcodes <- colnames(seu_anno)[seu_anno$cell_type == "Microglia" & seu_anno$orig.ident == "Mouse_Sham"]

seu_cf_mg   <- subset(seu_cf,   cells = colnames(seu_cf)[colnames(seu_cf) %in% strip_prefix(mg_cf_barcodes)])
seu_sham_mg <- subset(seu_sham, cells = colnames(seu_sham)[colnames(seu_sham) %in% strip_prefix(mg_sham_barcodes)])
seu_cf_mg$group   <- "Cystic_Fluid"
seu_sham_mg$group <- "Sham"

seu_mg <- merge(seu_cf_mg, y = seu_sham_mg, add.cell.ids = c("CF", "Sham"))
seu_mg <- JoinLayers(seu_mg)
rm(seu_cf_mg, seu_sham_mg); gc()

DefaultAssay(seu_mg) <- "RNA"
seu_mg <- NormalizeData(seu_mg, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
seu_mg <- FindVariableFeatures(seu_mg, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu_mg <- ScaleData(seu_mg, vars.to.regress = "percent.mt", verbose = FALSE)
seu_mg <- RunPCA(seu_mg, npcs = 30, verbose = FALSE)
seu_mg <- RunUMAP(seu_mg, dims = 1:20, verbose = FALSE)
seu_mg <- FindNeighbors(seu_mg, dims = 1:20, verbose = FALSE)
seu_mg <- FindClusters(seu_mg, resolution = 0.4, verbose = FALSE)

cat("\n=== CHECK 1: group table ===\n")
print(table(seu_mg$group))
cat("\n=== CHECK 2: cluster table ===\n")
print(table(seu_mg$seurat_clusters))

# Diagnostic UMAP
dir.create("results", showWarnings = FALSE, recursive = TRUE)
pdf("results/check_microglia_umap.pdf", width = 18, height = 12)
print(DimPlot(seu_mg, group.by = "group", pt.size = 0.5,
  cols = c("Cystic_Fluid" = "#E64B35", "Sham" = "#4DBBD5")) +
  ggtitle("CHECK: Microglia UMAP by Group"))
print(DimPlot(seu_mg, group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
  ggtitle("CHECK: Clusters"))
print(FeaturePlot(seu_mg, features = c("P2ry12","Tmem119","Cx3cr1","Sall1"),
  ncol = 4, pt.size = 0.4, order = TRUE) +
  plot_annotation(title = "CHECK: Homeostatic Markers"))
print(FeaturePlot(seu_mg, features = c("Apoe","Cd74","Trem2","Cst7"),
  ncol = 4, pt.size = 0.4, order = TRUE) +
  plot_annotation(title = "CHECK: Activated/DAM Markers"))
dev.off()

saveRDS(seu_mg, "data/processed/seu_microglia_merged.rds")

# ============================================================
# STEP 4-5: Build CDS & preprocess
# ============================================================
cat("\n=== STEP 4-5: Build CDS ===\n")
expr_mat <- GetAssayData(seu_mg, assay = "RNA", layer = "data")
gene_meta <- data.frame(gene_short_name = rownames(seu_mg), row.names = rownames(seu_mg))
cds <- new_cell_data_set(
  expression_data = expr_mat,
  cell_metadata   = seu_mg@meta.data,
  gene_metadata   = gene_meta
)
cat(sprintf("CDS: %d genes x %d cells\n", nrow(cds), ncol(cds)))

cds <- preprocess_cds(cds, num_dim = 20, method = "PCA", verbose = FALSE)
cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP", verbose = FALSE)

# ============================================================
# STEP 6: Cluster & learn graph (use_partition = TRUE)
# ============================================================
cat("\n=== STEP 6: Cluster & Graph ===\n")
cds <- cluster_cells(cds, resolution = 3e-4, verbose = TRUE)
cat(sprintf("Monocle3 clusters: %d\n", length(unique(monocle3::clusters(cds)))))

cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE, verbose = TRUE)
cat("Graph learned.\n")

# ============================================================
# STEP 7: Manual root selection (Sham homeostatic)
# ============================================================
cat("\n=== STEP 7: Manual Root ===\n")
root_markers <- c("P2ry12", "Tmem119", "Cx3cr1", "Sall1", "Hexb")
root_avail <- intersect(root_markers, rownames(cds))
cat(sprintf("Root markers: %d/5 available\n", length(root_avail)))

root_expr  <- colMeans(as.matrix(exprs(cds)[root_avail, , drop = FALSE]))
sham_cells <- colnames(cds)[colData(cds)$group == "Sham"]
sham_scores <- root_expr[sham_cells]
n_root <- max(10, round(length(sham_cells) * 0.05))
root_cells <- names(sort(sham_scores, decreasing = TRUE))[1:n_root]

cat(sprintf("Root cells: %d (all Sham, top 5%% homeostatic)\n", length(root_cells)))
stopifnot(all(root_cells %in% colnames(cds)))

# ============================================================
# STEP 8-9: Order cells + global normalization
# ============================================================
cat("\n=== STEP 8-9: Ordering + Normalization ===\n")
cds <- order_cells(cds, root_cells = root_cells)
pt <- pseudotime(cds)

cat(sprintf("Total: %d | Finite: %d | NA: %d | Inf: %d\n",
  length(pt), sum(is.finite(pt)), sum(is.na(pt)),
  sum(!is.finite(pt) & !is.na(pt))))
cat(sprintf("PT range: [%.3f, %.3f]\n", min(pt[is.finite(pt)]), max(pt[is.finite(pt)])))
cat("summary(pseudotime):\n"); print(summary(pt[is.finite(pt)]))

# Global normalization
fin <- is.finite(pt)
pt_min <- min(pt[fin]); pt_max <- max(pt[fin])
pt_norm <- rep(NA_real_, length(pt))
pt_norm[fin] <- (pt[fin] - pt_min) / (pt_max - pt_min)
colData(cds)$pseudotime_norm <- pt_norm
cat(sprintf("Norm PT range: [%.4f, %.4f]\n", min(pt_norm, na.rm = TRUE), max(pt_norm, na.rm = TRUE)))

# Per-group stats
for (g in c("Sham", "Cystic_Fluid")) {
  idx <- colData(cds)$group == g
  ptg <- pt_norm[idx]
  cat(sprintf("  %s: n_finite=%d, mean=%.4f, sd=%.4f\n", g, sum(is.finite(ptg)), mean(ptg, na.rm = TRUE), sd(ptg, na.rm = TRUE)))
}

# ============================================================
# FIGURES
# ============================================================

# ---- FIG ABC: UMAP + Trajectory ----
cat("\n=== Figures ===\n")
pdf("results/monocle3_fig_ABC_trajectory.pdf", width = 22, height = 7)
pA <- plot_cells(cds, color_cells_by = "group", label_cell_groups = FALSE, cell_size = 0.7) +
  scale_color_manual(values = c("Cystic_Fluid" = "#E64B35", "Sham" = "#4DBBD5")) +
  ggtitle("A: Microglia by Group")
pB <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2, cell_size = 0.7) +
  scale_color_viridis_c(option = "plasma") + ggtitle("B: Pseudotime on UMAP")
pC <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = TRUE,
  label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.6) +
  scale_color_viridis_c(option = "plasma") + ggtitle("C: Trajectory with Principal Graph")
print(pA + pB + pC + plot_layout(ncol = 3))
dev.off()

# ---- FIG D: Density ----
pdf("results/monocle3_fig_D_density.pdf", width = 8, height = 5)
pt_df <- data.frame(pt = pt_norm[fin], group = colData(cds)$group[fin])
wt <- wilcox.test(pt ~ group, data = pt_df)
pD <- ggplot(pt_df, aes(pt, fill = group)) +
  geom_density(alpha = 0.65) +
  scale_fill_manual(values = c("Cystic_Fluid" = "#E64B35", "Sham" = "#4DBBD5")) +
  labs(title = "D: Pseudotime Distribution (Sham vs CF)",
       subtitle = sprintf("Wilcoxon p = %.2e | Sham n=%d, CF n=%d",
         wt$p.value, sum(pt_df$group == "Sham"), sum(pt_df$group == "Cystic_Fluid")),
       x = "Normalized Pseudotime", y = "Density") +
  theme_minimal(base_size = 12)
print(pD)
dev.off()

# ---- FIG E: Homeostatic genes along pseudotime ----
pdf("results/monocle3_fig_E_homeostatic.pdf", width = 16, height = 4)
homeo <- intersect(c("P2ry12", "Tmem119", "Cx3cr1", "Sall1"), rownames(cds))
pE <- lapply(homeo, function(g) {
  plot_cells(cds, genes = g, label_cell_groups = FALSE,
    show_trajectory_graph = FALSE, cell_size = 0.5) +
    ggtitle(paste0("E: ", g)) + theme(legend.position = "none")
})
print(wrap_plots(pE, ncol = 4))
dev.off()

# ---- FIG F: Activated/DAM/antigen-presentation genes ----
pdf("results/monocle3_fig_F_activated.pdf", width = 18, height = 9)
dam <- intersect(c("Apoe", "Lpl", "Trem2", "Cst7", "Itgax", "Cd74", "H2-Aa", "H2-Ab1"), rownames(cds))
pF <- lapply(dam, function(g) {
  plot_cells(cds, genes = g, label_cell_groups = FALSE,
    show_trajectory_graph = FALSE, cell_size = 0.5) +
    ggtitle(paste0("F: ", g)) + theme(legend.position = "none")
})
print(wrap_plots(pF, ncol = 4))
dev.off()

# ---- Save CDS first (before heatmap, in case it fails) ----
saveRDS(cds, "data/processed/microglia_monocle3_final_cds.rds")

# ---- FIG G: Gene module heatmap ----
tryCatch({
  pdf("results/monocle3_fig_G_heatmap.pdf", width = 14, height = 10)
  fin_idx <- which(fin)
  early_cells <- colnames(cds)[fin_idx][pt_norm[fin_idx] < 0.25]
  late_cells  <- colnames(cds)[fin_idx][pt_norm[fin_idx] > 0.75]
  cat(sprintf("Heatmap: Early=%d cells, Late=%d cells\n", length(early_cells), length(late_cells)))

  if (length(early_cells) >= 10 && length(late_cells) >= 10) {
    early_mean <- rowMeans(as.matrix(exprs(cds)[, early_cells, drop = FALSE]))
    late_mean  <- rowMeans(as.matrix(exprs(cds)[, late_cells,  drop = FALSE]))
    fc <- late_mean - early_mean
    de_genes <- names(sort(abs(fc), decreasing = TRUE))[1:min(50, nrow(cds))]
    heatmap_cells <- c(early_cells, late_cells)
    hm_mat <- as.matrix(exprs(cds)[de_genes, heatmap_cells])
    hm_mat <- t(scale(t(hm_mat)))
    hm_mat[is.na(hm_mat)] <- 0
    hm_mat[hm_mat > 2] <- 2
    hm_mat[hm_mat < -2] <- -2

    col_ann <- data.frame(
      Group = as.character(colData(cds)[heatmap_cells, "group"]),
      Pseudotime = as.numeric(pt_norm[heatmap_cells]),
      row.names = heatmap_cells, stringsAsFactors = FALSE
    )

    library(pheatmap)
    pheatmap(hm_mat, cluster_cols = FALSE, cluster_rows = TRUE,
      show_colnames = FALSE, show_rownames = TRUE, fontsize_row = 6,
      annotation_col = col_ann,
      annotation_colors = list(
        Group = c(Cystic_Fluid = "#E64B35", Sham = "#4DBBD5"),
        Pseudotime = viridis::viridis(100)),
      main = "G: Top 50 DE Genes along Pseudotime (Early vs Late)",
      silent = TRUE)
  }
  dev.off()
  cat("Heatmap done.\n")
}, error = function(e) {
  cat("Heatmap failed:", conditionMessage(e), "\n")
  dev.off()
})

cat("\n========================================\n")
cat("ALL DONE. Output files:\n")
cat("  check_microglia_umap.pdf — diagnostic UMAPs\n")
cat("  monole3_fig_ABC_trajectory.pdf\n")
cat("  monole3_fig_D_density.pdf\n")
cat("  monole3_fig_E_homeostatic.pdf\n")
cat("  monole3_fig_F_activated.pdf\n")
cat("  monole3_fig_G_heatmap.pdf\n")
cat("========================================\n")
