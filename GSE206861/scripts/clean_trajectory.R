###############################################################################
# CLEAN trajectory — 从头构建, 只保留主分区, 单条干净轨迹
###############################################################################
library(Seurat)
library(monocle3)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(42)

# ============================================================
# PART 1: Build microglia Seurat (from raw QC data)
# ============================================================
cat("=== PART 1: Build microglia Seurat ===\n")
seu_cf <- readRDS("data/processed/seu_cf_qc.rds")
seu_sham <- readRDS("data/processed/seu_sham_qc.rds")
seu_anno <- readRDS("data/processed/seu_mouse_annotated.rds")

strip_prefix <- function(x) sub("^(CF_|Sham_)", "", x)
mg_cf_barcodes   <- colnames(seu_anno)[seu_anno$cell_type == "Microglia" & seu_anno$orig.ident == "Mouse_CF"]
mg_sham_barcodes <- colnames(seu_anno)[seu_anno$cell_type == "Microglia" & seu_anno$orig.ident == "Mouse_Sham"]

seu_cf_mg   <- subset(seu_cf,   cells = colnames(seu_cf)[colnames(seu_cf) %in% strip_prefix(mg_cf_barcodes)])
seu_sham_mg <- subset(seu_sham, cells = colnames(seu_sham)[colnames(seu_sham) %in% strip_prefix(mg_sham_barcodes)])
seu_cf_mg$group   <- "CF"
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

cat(sprintf("Merged microglia: %d cells (CF=%d, Sham=%d)\n",
  ncol(seu_mg), sum(seu_mg$group=="CF"), sum(seu_mg$group=="Sham")))

# ============================================================
# PART 2: Build CDS
# ============================================================
cat("\n=== PART 2: Build CDS ===\n")
expr_mat <- GetAssayData(seu_mg, assay = "RNA", layer = "data")
gene_meta <- data.frame(gene_short_name = rownames(seu_mg), row.names = rownames(seu_mg))
cds <- new_cell_data_set(expression_data = expr_mat, cell_metadata = seu_mg@meta.data, gene_metadata = gene_meta)
cat(sprintf("CDS: %d genes x %d cells\n", nrow(cds), ncol(cds)))

# ============================================================
# PART 3: monocle3 — cluster + learn_graph with partitions
# ============================================================
cat("\n=== PART 3: monocle3 pipeline ===\n")
cds <- preprocess_cds(cds, num_dim = 20, method = "PCA", verbose = FALSE)
cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP", verbose = FALSE)
cds <- cluster_cells(cds, resolution = 3e-4, verbose = TRUE)
cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE, verbose = FALSE)

partitions <- monocle3::partitions(cds)
cat("Partitions found:", length(unique(partitions)), "\n")
cat("Partition sizes:", paste(table(partitions), collapse=", "), "\n")

# ============================================================
# PART 4: Select main partition (must contain both groups + homeostatic signal)
# ============================================================
cat("\n=== PART 4: Select main partition ===\n")
root_markers <- c("P2ry12","Tmem119","Cx3cr1","Sall1","Hexb")
root_avail <- intersect(root_markers, rownames(cds))
root_expr <- colMeans(as.matrix(exprs(cds)[root_avail, , drop = FALSE]))

best_p <- NULL; best_score <- -Inf
for(p in sort(unique(partitions))) {
  cells_p <- colnames(cds)[partitions == p]
  grp <- colData(cds)[cells_p, "group"]
  n_sham <- as.integer(sum(grp == "Sham"))
  n_cf   <- as.integer(sum(grp == "CF"))
  sham_in_p <- cells_p[grp == "Sham"]
  homeo_score <- if(length(sham_in_p) > 5) mean(root_expr[sham_in_p]) else 0
  cat(sprintf("  P%s: %d cells (Sham=%d, CF=%d) homeo=%.3f\n",
    as.character(p), length(cells_p), n_sham, n_cf, homeo_score))
  # Must have both groups and high homeostatic score
  if(n_sham >= 5 && n_cf >= 5 && homeo_score > 0.3 && homeo_score > best_score) {
    best_p <- p; best_score <- homeo_score
  }
}
if(is.null(best_p)) best_p <- 1  # fallback to largest
cat(sprintf("Selected partition: %s (homeo_score=%.3f)\n", as.character(best_p), best_score))

# ============================================================
# PART 5: Subset & re-learn as single connected graph
# ============================================================
cat("\n=== PART 5: Re-learn clean graph ===\n")
cds <- cds[, partitions == best_p]
cat(sprintf("Subset to %d cells\n", ncol(cds)))

# CRITICAL: re-cluster after subsetting to refresh monocle3 internal state
cds <- cluster_cells(cds, resolution = 3e-4, verbose = TRUE)
cat(sprintf("Re-clustered: %d clusters\n", length(unique(monocle3::clusters(cds)))))

# Now safe to use use_partition=FALSE — all cells in same partition
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE, verbose = TRUE)

# ============================================================
# PART 6: Manual root (Sham homeostatic)
# ============================================================
cat("\n=== PART 6: Root selection ===\n")
root_expr_sub <- colMeans(as.matrix(exprs(cds)[root_avail, , drop = FALSE]))
sham_cells <- colnames(cds)[colData(cds)$group == "Sham"]
sham_scores <- root_expr_sub[sham_cells]
n_root <- max(10, round(length(sham_cells) * 0.05))
root_cells <- names(sort(sham_scores, decreasing = TRUE))[1:n_root]
cat(sprintf("Root: %d cells (all Sham, score=%.3f)\n", length(root_cells), mean(root_expr_sub[root_cells])))

# ============================================================
# PART 7: Order + normalize
# ============================================================
cds <- order_cells(cds, root_cells = root_cells)
pt <- pseudotime(cds)
cat(sprintf("\nPseudotime checks:\n  Total=%d, Finite=%d, NA=%d, Inf=%d\n",
  length(pt), sum(is.finite(pt)), sum(is.na(pt)), sum(!is.finite(pt) & !is.na(pt))))
cat(sprintf("  Range: [%.3f, %.3f]\n", min(pt[is.finite(pt)]), max(pt[is.finite(pt)])))

fin <- is.finite(pt)
pt_norm <- rep(NA_real_, length(pt))
pt_norm[fin] <- (pt[fin] - min(pt[fin])) / (max(pt[fin]) - min(pt[fin]))
colData(cds)$pseudotime_norm <- pt_norm

for(g in c("Sham", "CF")) {
  idx <- colData(cds)$group == g
  cat(sprintf("  %s: n=%d, mean_pt=%.4f, sd=%.4f\n", g, sum(idx), mean(pt_norm[idx], na.rm=TRUE), sd(pt_norm[idx], na.rm=TRUE)))
}

# ============================================================
# PART 8: FIGURES
# ============================================================
cat("\n=== PART 8: Figures ===\n")

# --- A+B+C: UMAP + Trajectory ---
pdf("results/monocle3_CLEAN_trajectory.pdf", width = 22, height = 7)
pA <- plot_cells(cds, color_cells_by = "group", label_cell_groups = FALSE, cell_size = 0.8) +
  scale_color_manual(values = c("CF" = "#E64B35", "Sham" = "#4DBBD5")) +
  ggtitle("A: Microglia by Group")
pB <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2, cell_size = 0.8) +
  scale_color_viridis_c(option = "plasma") + ggtitle("B: Pseudotime")
pC <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = TRUE,
  label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.6) +
  scale_color_viridis_c(option = "plasma") + ggtitle("C: Trajectory Graph")
print(pA + pB + pC + plot_layout(ncol = 3))
dev.off()

# --- D: Density ---
pdf("results/monocle3_CLEAN_density.pdf", width = 8, height = 5)
pt_df <- data.frame(pt = pt_norm[fin], group = colData(cds)$group[fin])
wt <- wilcox.test(pt ~ group, data = pt_df)
print(ggplot(pt_df, aes(pt, fill = group)) +
  geom_density(alpha = 0.65) +
  scale_fill_manual(values = c("CF" = "#E64B35", "Sham" = "#4DBBD5")) +
  labs(title = "D: Pseudotime Distribution",
       subtitle = sprintf("Wilcoxon p = %.2e | Sham=%.3f, CF=%.3f",
         wt$p.value, mean(pt_df$pt[pt_df$group=="Sham"]), mean(pt_df$pt[pt_df$group=="CF"])),
       x = "Normalized Pseudotime", y = "Density") +
  theme_minimal(base_size = 12))
dev.off()

# --- E+F: All marker genes ---
pdf("results/monocle3_CLEAN_genes.pdf", width = 18, height = 14)
homeo <- intersect(c("P2ry12","Tmem119","Cx3cr1","Sall1"), rownames(cds))
dam   <- intersect(c("Apoe","Lpl","Trem2","Cst7","Itgax","Cd74","H2-Aa","H2-Ab1"), rownames(cds))
all_g <- c(homeo, dam)
p_list <- lapply(all_g, function(g) {
  plot_cells(cds, genes = g, label_cell_groups = FALSE,
    show_trajectory_graph = FALSE, cell_size = 0.5) +
    ggtitle(g) + theme(legend.position = "none")
})
print(wrap_plots(p_list, ncol = 4))
dev.off()

# --- G: Heatmap (DE early vs late) ---
pdf("results/monocle3_CLEAN_heatmap.pdf", width = 14, height = 10)
fin_idx <- which(fin)
early <- colnames(cds)[fin_idx][pt_norm[fin_idx] < 0.25]
late  <- colnames(cds)[fin_idx][pt_norm[fin_idx] > 0.75]
cat(sprintf("Heatmap: early=%d, late=%d\n", length(early), length(late)))

if(length(early)>=10 && length(late)>=10) {
  fc_vec <- rowMeans(as.matrix(exprs(cds)[, late, drop=FALSE])) -
            rowMeans(as.matrix(exprs(cds)[, early, drop=FALSE]))
  de_genes <- names(sort(abs(fc_vec), decreasing = TRUE))[1:50]
  hm_cells <- c(early, late)
  hm_mat <- as.matrix(exprs(cds)[de_genes, hm_cells])
  hm_mat <- t(scale(t(hm_mat)))
  hm_mat[is.na(hm_mat)] <- 0
  hm_mat[hm_mat > 2] <- 2; hm_mat[hm_mat < -2] <- -2

  col_ann <- data.frame(
    Group = colData(cds)[hm_cells, "group"],
    row.names = hm_cells, stringsAsFactors = FALSE
  )
  col_ann$Group <- factor(col_ann$Group, levels = c("CF", "Sham"))

  library(pheatmap)
  pheatmap(hm_mat, cluster_cols = FALSE, cluster_rows = TRUE,
    show_colnames = FALSE, show_rownames = TRUE, fontsize_row = 6,
    annotation_col = col_ann,
    annotation_colors = list(Group = c(CF = "#E64B35", Sham = "#4DBBD5")),
    gaps_col = length(early),
    main = "G: Top 50 DE Genes (Early vs Late Pseudotime)",
    border_color = NA)
}
dev.off()

# ============================================================
# Save
# ============================================================
saveRDS(cds, "data/processed/microglia_monocle3_clean.rds")

cat("\n===========================================\n")
cat("CLEAN TRAJECTORY COMPLETE!\n")
cat("Output:\n")
cat("  monole3_CLEAN_trajectory.pdf (A+B+C)\n")
cat("  monole3_CLEAN_density.pdf   (D)\n")
cat("  monole3_CLEAN_genes.pdf     (E+F)\n")
cat("  monole3_CLEAN_heatmap.pdf   (G)\n")
cat("===========================================\n")
