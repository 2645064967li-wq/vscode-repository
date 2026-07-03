###############################################################################
# GSE206861 单细胞分析 — 标准化、整合、降维、聚类
# 小鼠样本 (CF vs Sham) 整合；人类样本单独分析
###############################################################################

library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(here)
library(future)

# ============================================================================
# 0. 环境设置
# ============================================================================
setwd(here("GSE206861"))

# 并行计算设置
plan("multisession", workers = max(1, parallel::detectCores() - 2))
options(future.globals.maxSize = 8 * 1024^3)  # 8GB
set.seed(42)

# 读取 QC 后数据
seu_cf   <- readRDS(here("GSE206861/data/processed/seu_cf_qc.rds"))
seu_sham <- readRDS(here("GSE206861/data/processed/seu_sham_qc.rds"))
seu_acp  <- readRDS(here("GSE206861/data/processed/seu_acp_qc.rds"))

# ============================================================================
# Part A: 小鼠分析 — 整合 + 聚类
# ============================================================================

cat("\n========== Part A: 小鼠整合分析 ==========\n")

# --- A1. 合并小鼠样本 ---
seu_mouse <- merge(seu_cf, y = seu_sham, add.cell.ids = c("CF", "Sham"))
# Clean up unused levels and metadata
seu_mouse$orig.ident <- factor(seu_mouse$orig.ident, levels = c("Mouse_CF", "Mouse_Sham"))

cat("Merged mouse object:", ncol(seu_mouse), "cells,", nrow(seu_mouse), "genes\n")

# --- A2. SCTransform v2 标准化（替代 LogNormalize + ScaleData + FindVariableFeatures）---
# 注意：SCTransform 本身包含了标准化、高变基因选择、方差稳定化
cat("\n>>> SCTransform v2 标准化...\n")

seu_mouse <- SCTransform(
  seu_mouse,
  vst.flavor = "v2",
  vars.to.regress = c("percent.mt"),
  verbose = TRUE
)

# --- A3. Run PCA ---
cat("\n>>> PCA...\n")
seu_mouse <- RunPCA(seu_mouse, npcs = 50, verbose = FALSE)

# 查看 PCA 维度选择
pdf(here("GSE206861/results/mouse_elbow.pdf"), width = 8, height = 5)
ElbowPlot(seu_mouse, ndims = 50) + ggtitle("Mouse: PCA Elbow Plot")
dev.off()

# --- A4. Harmony 整合（校正样本间批次效应）---
cat("\n>>> Harmony 整合...\n")
seu_mouse <- RunHarmony(
  seu_mouse,
  group.by.vars = "orig.ident",
  reduction = "pca",
  assay.use = "SCT",
  dims.use = 1:30,
  theta = 2,
  sigma = 0.1,
  max.iter.harmony = 20
)

# --- A5. UMAP 降维 ---
cat("\n>>> UMAP...\n")
seu_mouse <- RunUMAP(seu_mouse, reduction = "harmony", dims = 1:30, verbose = FALSE)

# --- A6. 构建 SNN 图 + 聚类 ---
cat("\n>>> 聚类...\n")
seu_mouse <- FindNeighbors(seu_mouse, reduction = "harmony", dims = 1:30)
seu_mouse <- FindClusters(seu_mouse, resolution = seq(0.2, 1.2, by = 0.2))

# 选择最佳分辨率 (使用 0.6 作为默认)
Idents(seu_mouse) <- "SCT_snn_res.0.6"
seu_mouse$seurat_clusters <- seu_mouse$SCT_snn_res.0.6

# --- A7. 聚类树状图（评估不同分辨率）---
pdf(here("GSE206861/results/mouse_clustree.pdf"), width = 12, height = 8)
if (requireNamespace("clustree", quietly = TRUE)) {
  library(clustree)
  clustree(seu_mouse, prefix = "SCT_snn_res.")
} else {
  cat("  (安装 clustree 可绘制聚类树: install.packages('clustree'))\n")
}
dev.off()

# --- A8. UMAP 可视化 ---
p_cluster <- DimPlot(seu_mouse, group.by = "seurat_clusters", label = TRUE, repel = TRUE,
                     pt.size = 0.3) + ggtitle("Mouse: Clusters") + NoLegend()
p_sample  <- DimPlot(seu_mouse, group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("Mouse: Sample (CF vs Sham)") + scale_color_manual(values = c("Mouse_CF" = "#E64B35", "Mouse_Sham" = "#4DBBD5"))

pdf(here("GSE206861/results/mouse_umap.pdf"), width = 16, height = 7)
p_cluster + p_sample
dev.off()

cat("\n  Mouse clustering: ", length(unique(seu_mouse$seurat_clusters)), "clusters\n")

# --- A9. 保存 ---
saveRDS(seu_mouse, here("GSE206861/data/processed/seu_mouse_integrated.rds"))
rm(seu_cf, seu_sham)
gc()

# ============================================================================
# Part B: 人类样本 — 单独降维聚类
# ============================================================================

cat("\n========== Part B: 人类 ACP 分析 ==========\n")

# --- B1. SCTransform ---
cat("\n>>> SCTransform v2...\n")
seu_acp <- SCTransform(
  seu_acp,
  vst.flavor = "v2",
  vars.to.regress = c("percent.mt"),
  verbose = TRUE
)

# --- B2. PCA + UMAP + 聚类 ---
cat("\n>>> PCA + UMAP + Clustering...\n")
seu_acp <- RunPCA(seu_acp, npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = seq(0.2, 1.0, by = 0.2))

Idents(seu_acp) <- "SCT_snn_res.0.6"
seu_acp$seurat_clusters <- seu_acp$SCT_snn_res.0.6

# --- B3. UMAP ---
pdf(here("GSE206861/results/human_umap.pdf"), width = 8, height = 7)
DimPlot(seu_acp, group.by = "seurat_clusters", label = TRUE, repel = TRUE,
        pt.size = 0.5) + ggtitle("Human ACP: Clusters")
dev.off()

cat("  Human ACP clustering:", length(unique(seu_acp$seurat_clusters)), "clusters\n")

# --- B4. 保存 ---
saveRDS(seu_acp, here("GSE206861/data/processed/seu_acp_integrated.rds"))

# ============================================================================
# 完成
# ============================================================================
cat("\n✓ 标准化、整合、降维、聚类完成！\n")
cat("  保存: data/processed/seu_mouse_integrated.rds\n")
cat("  保存: data/processed/seu_acp_integrated.rds\n")
cat("  接下来运行: 04_cell_annotation.R\n")

plan("sequential")  # 清理并行设置
