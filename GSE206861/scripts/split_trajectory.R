###############################################################################
# 两组独立拟时序 — CF vs Sham 各自分化轨迹, 并排比较
###############################################################################
library(Seurat)
library(monocle3)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleCellExperiment)
set.seed(42)

# ============================================================
# PART 1: Build separate Seurat objects for CF and Sham microglia
# ============================================================
cat("=== PART 1: Extract microglia ===\n")

seu_cf   <- readRDS("data/processed/seu_cf_qc.rds")
seu_sham <- readRDS("data/processed/seu_sham_qc.rds")
seu_anno <- readRDS("data/processed/seu_mouse_annotated.rds")

strip_prefix <- function(x) sub("^(CF_|Sham_)", "", x)

mg_cf_barcodes   <- colnames(seu_anno)[seu_anno$cell_type == "Microglia" & seu_anno$orig.ident == "Mouse_CF"]
mg_sham_barcodes <- colnames(seu_anno)[seu_anno$cell_type == "Microglia" & seu_anno$orig.ident == "Mouse_Sham"]

seu_cf_mg   <- subset(seu_cf,   cells = colnames(seu_cf)[colnames(seu_cf) %in% strip_prefix(mg_cf_barcodes)])
seu_sham_mg <- subset(seu_sham, cells = colnames(seu_sham)[colnames(seu_sham) %in% strip_prefix(mg_sham_barcodes)])

cat(sprintf("CF: %d cells, Sham: %d cells\n", ncol(seu_cf_mg), ncol(seu_sham_mg)))

# ============================================================
# PART 2: Process each Seurat object independently
# ============================================================
process_microglia <- function(seu_obj, label) {
  cat(sprintf("\n--- Processing %s ---\n", label))
  DefaultAssay(seu_obj) <- "RNA"
  seu_obj <- NormalizeData(seu_obj, method="LogNormalize", scale.factor=10000, verbose=FALSE)
  seu_obj <- FindVariableFeatures(seu_obj, method="vst", nfeatures=2000, verbose=FALSE)
  seu_obj <- ScaleData(seu_obj, vars.to.regress="percent.mt", verbose=FALSE)
  seu_obj <- RunPCA(seu_obj, npcs=20, verbose=FALSE)
  seu_obj <- RunUMAP(seu_obj, dims=1:15, verbose=FALSE)
  seu_obj <- FindNeighbors(seu_obj, dims=1:15, verbose=FALSE)
  seu_obj <- FindClusters(seu_obj, resolution=0.3, verbose=FALSE)
  cat(sprintf("  Clusters: %d\n", length(unique(seu_obj$seurat_clusters))))
  return(seu_obj)
}

seu_cf_mg   <- process_microglia(seu_cf_mg, "CF")
seu_sham_mg <- process_microglia(seu_sham_mg, "Sham")

# ============================================================
# PART 3: Build CDS and run monocle3
# ============================================================
run_monocle3 <- function(seu_obj, label, root_markers=c("P2ry12","Tmem119","Cx3cr1","Sall1","Hexb")) {
  cat(sprintf("\n=== monocle3: %s ===\n", label))

  expr_mat <- GetAssayData(seu_obj, assay="RNA", layer="data")
  gene_meta <- data.frame(gene_short_name=rownames(seu_obj), row.names=rownames(seu_obj))
  cds <- new_cell_data_set(expression_data=expr_mat, cell_metadata=seu_obj@meta.data, gene_metadata=gene_meta)
  cat(sprintf("  CDS: %d genes x %d cells\n", nrow(cds), ncol(cds)))

  cds <- preprocess_cds(cds, num_dim=20, method="PCA", verbose=FALSE)
  cds <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="UMAP", verbose=FALSE)
  cds <- cluster_cells(cds, resolution=3e-4, verbose=FALSE)

  # Check partitions
  cds <- learn_graph(cds, use_partition=TRUE, close_loop=FALSE, verbose=FALSE)
  partitions <- monocle3::partitions(cds)
  cat(sprintf("  Partitions: %d (sizes: %s)\n", length(unique(partitions)),
    paste(table(partitions), collapse=", ")))

  # Select main partition
  root_avail <- intersect(root_markers, rownames(cds))
  root_expr <- colMeans(as.matrix(exprs(cds)[root_avail, , drop=FALSE]))

  best_p <- NULL; best_score <- -Inf
  for(p in sort(unique(partitions))) {
    cells_p <- colnames(cds)[partitions == p]
    homeo_p <- mean(root_expr[cells_p])
    n_p <- length(cells_p)
    if(n_p > 20 && homeo_p > best_score) { best_p <- p; best_score <- homeo_p }
  }
  if(is.null(best_p)) best_p <- names(which.max(table(partitions)))
  cat(sprintf("  Selected partition %s (%d cells, homeo_score=%.3f)\n",
    as.character(best_p), sum(partitions==best_p), best_score))

  # Subset to main partition & re-learn clean graph
  cds <- cds[, partitions == best_p]
  cds <- cluster_cells(cds, resolution=3e-4, verbose=FALSE)
  cds <- learn_graph(cds, use_partition=FALSE, close_loop=FALSE, verbose=FALSE)

  # Root = highest homeostatic marker expression
  root_expr_sub <- colMeans(as.matrix(exprs(cds)[root_avail, , drop=FALSE]))
  n_root <- max(5, round(ncol(cds) * 0.03))
  root_cells <- names(sort(root_expr_sub, decreasing=TRUE))[1:n_root]
  cat(sprintf("  Root: %d cells, score=%.3f\n", length(root_cells), mean(root_expr_sub[root_cells])))

  cds <- order_cells(cds, root_cells=root_cells)
  pt <- pseudotime(cds)
  cat(sprintf("  Pseudotime: %d/%d finite, range [%.2f, %.2f]\n",
    sum(is.finite(pt)), length(pt),
    min(pt[is.finite(pt)]), max(pt[is.finite(pt)])))

  # Global normalization (within this group)
  fin <- is.finite(pt)
  pt_norm <- rep(NA_real_, length(pt))
  pt_norm[fin] <- (pt[fin] - min(pt[fin])) / (max(pt[fin]) - min(pt[fin]))
  colData(cds)$pseudotime_norm <- pt_norm
  cat(sprintf("  Normalized PT: mean=%.4f, sd=%.4f\n", mean(pt_norm,na.rm=TRUE), sd(pt_norm,na.rm=TRUE)))

  # Rename UMAP dims
  umap_coords <- reducedDim(cds, "UMAP")
  colnames(umap_coords) <- c("Component 1", "Component 2")
  reducedDim(cds, "UMAP") <- umap_coords

  return(cds)
}

cds_cf   <- run_monocle3(seu_cf_mg, "CF")
cds_sham <- run_monocle3(seu_sham_mg, "Sham")

# ============================================================
# PART 4: Side-by-side comparison figures
# ============================================================

# --- FIG 1: Trajectory comparison ---
pdf("results/monocle3_split_trajectory.pdf", width=22, height=14)

# Row 1: Trajectories
pCF_traj <- plot_cells(cds_cf, color_cells_by="pseudotime", label_cell_groups=TRUE,
  label_leaves=TRUE, label_branch_points=TRUE, graph_label_size=3, cell_size=0.7) +
  scale_color_viridis_c(option="plasma") +
  ggtitle("CF (Cystic Fluid): Microglia Trajectory",
          subtitle=sprintf("%d cells | PT range [%.2f, %.2f]",
            ncol(cds_cf), min(pseudotime(cds_cf),na.rm=TRUE), max(pseudotime(cds_cf),na.rm=TRUE)))

pSham_traj <- plot_cells(cds_sham, color_cells_by="pseudotime", label_cell_groups=TRUE,
  label_leaves=TRUE, label_branch_points=TRUE, graph_label_size=3, cell_size=0.7) +
  scale_color_viridis_c(option="plasma") +
  ggtitle("Sham (PBS): Microglia Trajectory",
          subtitle=sprintf("%d cells | PT range [%.2f, %.2f]",
            ncol(cds_sham), min(pseudotime(cds_sham),na.rm=TRUE), max(pseudotime(cds_sham),na.rm=TRUE)))

# Row 2: Pseudotime density
pt_cf <- pseudotime(cds_cf); pt_cf <- pt_cf[is.finite(pt_cf)]
pt_sham <- pseudotime(cds_sham); pt_sham <- pt_sham[is.finite(pt_sham)]

# Normalize BOTH together for fair comparison
# (min across both = 0, max across both = 1)
all_pt <- c(pt_cf, pt_sham)
global_min <- min(all_pt); global_max <- max(all_pt)
cf_norm <- (pt_cf - global_min) / (global_max - global_min)
sham_norm <- (pt_sham - global_min) / (global_max - global_min)

pt_df <- rbind(
  data.frame(pt=cf_norm, Group="CF"),
  data.frame(pt=sham_norm, Group="Sham")
)

pDensity <- ggplot(pt_df, aes(pt, fill=Group)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values=c("CF"="#E64B35","Sham"="#4DBBD5")) +
  labs(title="Pseudotime Distribution Comparison",
       subtitle=sprintf("CF mean=%.3f | Sham mean=%.3f | Wilcoxon p=%.2e",
         mean(cf_norm), mean(sham_norm),
         wilcox.test(pt~Group, data=pt_df)$p.value),
       x="Normalized Pseudotime (0=root/homeostatic)", y="Density") +
  theme_minimal(base_size=12)

# Row 3: Trajectory complexity metrics
cf_n_branches <- length(unique(monocle3::clusters(cds_cf)))
sham_n_branches <- length(unique(monocle3::clusters(cds_sham)))

pSummary <- ggplot() +
  annotate("text", x=1, y=3, label=sprintf("CF:  %d cells, %d branches, PT range [%.1f, %.1f]",
    ncol(cds_cf), cf_n_branches, min(pt_cf), max(pt_cf)), size=5, hjust=0, color="#E64B35") +
  annotate("text", x=1, y=2, label=sprintf("Sham: %d cells, %d branches, PT range [%.1f, %.1f]",
    ncol(cds_sham), sham_n_branches, min(pt_sham), max(pt_sham)), size=5, hjust=0, color="#4DBBD5") +
  annotate("text", x=1, y=1, label="Both rooted at P2ry12/Tmem119/Cx3cr1/Sall1/Hexb high cells",
    size=4, hjust=0, color="grey40") +
  xlim(1, 10) + ylim(0, 4) + theme_void() +
  ggtitle("Trajectory Summary")

print(pCF_traj + pSham_traj + plot_layout(ncol=2))
print(pDensity + pSummary + plot_layout(ncol=2))
dev.off()

# --- FIG 2: Gene expression comparison ---
pdf("results/monocle3_split_genes.pdf", width=20, height=16)

cf_genes <- intersect(c("P2ry12","Tmem119","Cx3cr1","Sall1",
                        "Apoe","Lpl","Trem2","Cst7","Itgax","Cd74","H2-Aa","H2-Ab1",
                        "C1qa","C1qb","C1qc"), rownames(cds_cf))
sham_genes <- intersect(c("P2ry12","Tmem119","Cx3cr1","Sall1",
                          "Apoe","Lpl","Trem2","Cst7","Itgax","Cd74","H2-Aa","H2-Ab1",
                          "C1qa","C1qb","C1qc"), rownames(cds_sham))
common_genes <- intersect(cf_genes, sham_genes)

for(g in common_genes) {
  p_g_cf <- plot_cells(cds_cf, genes=g, label_cell_groups=FALSE,
    show_trajectory_graph=FALSE, cell_size=0.5) +
    ggtitle(sprintf("CF: %s", g)) + theme(legend.position="none")
  p_g_sham <- plot_cells(cds_sham, genes=g, label_cell_groups=FALSE,
    show_trajectory_graph=FALSE, cell_size=0.5) +
    ggtitle(sprintf("Sham: %s", g)) + theme(legend.position="none")
  print(p_g_cf + p_g_sham + plot_layout(ncol=2))
}
dev.off()

# --- FIG 3: Cd74 head-to-head spotlight ---
pdf("results/monocle3_split_Cd74_spotlight.pdf", width=14, height=6)
p_cf_cd74 <- plot_cells(cds_cf, genes="Cd74", label_cell_groups=FALSE,
  show_trajectory_graph=TRUE, cell_size=0.7) +
  scale_color_viridis_c(option="inferno") +
  ggtitle("CF: Cd74 (CD74-APP axis)")
p_sham_cd74 <- plot_cells(cds_sham, genes="Cd74", label_cell_groups=FALSE,
  show_trajectory_graph=TRUE, cell_size=0.7) +
  scale_color_viridis_c(option="inferno") +
  ggtitle("Sham: Cd74")
print(p_cf_cd74 + p_sham_cd74 + plot_layout(ncol=2))
dev.off()

# Save
saveRDS(cds_cf, "data/processed/microglia_cf_monocle3.rds")
saveRDS(cds_sham, "data/processed/microglia_sham_monocle3.rds")

cat("\n===========================================\n")
cat("SPLIT TRAJECTORY COMPLETE!\n")
cat("  monole3_split_trajectory.pdf\n")
cat("  monole3_split_genes.pdf\n")
cat("  monole3_split_Cd74_spotlight.pdf\n")
cat("===========================================\n")
