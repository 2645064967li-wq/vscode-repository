###############################################################################
# GSE206861 单细胞分析 — 数据读入与质控
# 文献: ACP cyst fluid triggers microglia activation (PMID: 35525962)
# 3 个样本: Mouse_Cystic-fluid, Mouse_Sham, Human_ACP
###############################################################################

# ============================================================================
# 0. 环境准备
# ============================================================================
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(here)

# 设置工作目录
setwd(here("GSE206861"))
data_dir <- here("GSE206861/data")

# ============================================================================
# 1. 数据读入 — 10x Genomics 三元格式
# ============================================================================

cat("\n========== 读入数据 ==========\n")

# --- GSM6265811: Mouse Cystic-fluid (下丘脑注射 ACP 囊液) ---
cat("\n>>> 读入 Mouse_Cystic-fluid...\n")
seu_cf <- Read10X(data.dir = file.path(data_dir, "GSM6265811_Mouse_Cystic-fluid"))
seu_cf <- CreateSeuratObject(
  counts  = seu_cf,
  project = "Mouse_CF",
  assay   = "RNA",
  min.cells = 3,
  min.features = 200
)
seu_cf$sample <- "Cystic-fluid"
seu_cf$species <- "Mouse"
cat("  Cells:", ncol(seu_cf), "| Genes:", nrow(seu_cf), "\n")

# --- GSM6265812: Mouse Sham (PBS 对照) ---
cat("\n>>> 读入 Mouse_Sham...\n")
seu_sham <- Read10X(data.dir = file.path(data_dir, "GSM6265812_Mouse_Sham"))
seu_sham <- CreateSeuratObject(
  counts  = seu_sham,
  project = "Mouse_Sham",
  assay   = "RNA",
  min.cells = 3,
  min.features = 200
)
seu_sham$sample <- "Sham"
seu_sham$species <- "Mouse"
cat("  Cells:", ncol(seu_sham), "| Genes:", nrow(seu_sham), "\n")

# --- GSM6265813: Human ACP tumor tissue ---
cat("\n>>> 读入 Human_ACP...\n")
seu_acp <- Read10X(data.dir = file.path(data_dir, "GSM6265813_Human_ACP"))
seu_acp <- CreateSeuratObject(
  counts  = seu_acp,
  project = "Human_ACP",
  assay   = "RNA",
  min.cells = 3,
  min.features = 200
)
seu_acp$sample <- "Human_ACP"
seu_acp$species <- "Human"
cat("  Cells:", ncol(seu_acp), "| Genes:", nrow(seu_acp), "\n")

# ============================================================================
# 2. 计算 QC 指标
# ============================================================================

cat("\n========== 计算 QC 指标 ==========\n")

# --- 通用 QC 指标函数 ---
add_qc_metrics <- function(seu, species = c("Mouse", "Human")) {
  species <- match.arg(species)

  # 线粒体基因
  if (species == "Mouse") {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
    seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rp[sl]")
    seu[["percent.hb"]] <- PercentageFeatureSet(seu, pattern = "^Hb[^(p)]")
  } else {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
    seu[["percent.hb"]] <- PercentageFeatureSet(seu, pattern = "^HB[^(P)]")
  }

  # log10(genes per UMI) — 复杂度
  seu[["log10GenesPerUMI"]] <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)

  return(seu)
}

seu_cf   <- add_qc_metrics(seu_cf,   "Mouse")
seu_sham <- add_qc_metrics(seu_sham, "Mouse")
seu_acp  <- add_qc_metrics(seu_acp,  "Human")

# ============================================================================
# 3. QC 可视化（过滤前）
# ============================================================================

cat("\n========== QC 可视化（过滤前）==========\n")

# --- 通用 Violin 图函数 ---
plot_qc_violin <- function(seu, title = "") {
  VlnPlot(seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    ncol = 4, pt.size = 0.1, group.by = "orig.ident"
  ) + plot_annotation(title = title)
}

# --- 散点图（nCount vs nFeature, 标注 %MT）---
plot_qc_scatter <- function(seu, title = "") {
  p1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        group.by = "orig.ident", pt.size = 0.5) +
    ggtitle(paste(title, ": Count vs Feature"))
  p2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        group.by = "orig.ident", pt.size = 0.5) +
    ggtitle(paste(title, ": Count vs %Mito"))
  p1 + p2
}

# 绘制
pdf(here("GSE206861/results/QC_prefilter_violin.pdf"), width = 16, height = 6)
plot_qc_violin(seu_cf, "Mouse Cystic-fluid (Pre-filter)")
plot_qc_violin(seu_sham, "Mouse Sham (Pre-filter)")
plot_qc_violin(seu_acp, "Human ACP (Pre-filter)")
dev.off()

pdf(here("GSE206861/results/QC_prefilter_scatter.pdf"), width = 12, height = 5)
plot_qc_scatter(seu_cf, "Mouse CF")
plot_qc_scatter(seu_sham, "Mouse Sham")
plot_qc_scatter(seu_acp, "Human ACP")
dev.off()

cat("  QC 图已保存至 results/QC_prefilter_*.pdf\n")

# ============================================================================
# 4. 设定 QC 过滤阈值
# ============================================================================

cat("\n========== QC 过滤 ==========\n")

# ----- 小鼠阈值 -----
qc_mouse <- list(
  nFeature_low   = 200,
  nFeature_high  = 6000,
  nCount_low     = 500,
  nCount_high    = 50000,
  percent.mt_max = 20,
  percent.ribo_max = 40
)

# ----- 人类阈值 -----
qc_human <- list(
  nFeature_low   = 200,
  nFeature_high  = 5000,
  nCount_low     = 500,
  nCount_high    = 40000,
  percent.mt_max = 20,
  percent.ribo_max = 40
)

# --- 应用过滤 ---
filter_cells <- function(seu, qc) {
  n_before <- ncol(seu)

  seu <- subset(seu,
    nFeature_RNA  > qc$nFeature_low  & nFeature_RNA < qc$nFeature_high &
    nCount_RNA    > qc$nCount_low    & nCount_RNA   < qc$nCount_high  &
    percent.mt    < qc$percent.mt_max &
    percent.ribo  < qc$percent.ribo_max
  )

  n_after <- ncol(seu)
  pct_removed <- round((1 - n_after / n_before) * 100, 1)
  cat("  Cells before:", n_before, "| after:", n_after, "| removed:", pct_removed, "%\n")

  return(seu)
}

cat("\n--- Mouse Cystic-fluid ---\n")
seu_cf <- filter_cells(seu_cf, qc_mouse)
cat("\n--- Mouse Sham ---\n")
seu_sham <- filter_cells(seu_sham, qc_mouse)
cat("\n--- Human ACP ---\n")
seu_acp <- filter_cells(seu_acp, qc_human)

# ============================================================================
# 5. QC 可视化（过滤后）
# ============================================================================

pdf(here("GSE206861/results/QC_postfilter_violin.pdf"), width = 16, height = 6)
plot_qc_violin(seu_cf, "Mouse Cystic-fluid (Post-filter)")
plot_qc_violin(seu_sham, "Mouse Sham (Post-filter)")
plot_qc_violin(seu_acp, "Human ACP (Post-filter)")
dev.off()

cat("\n  QC 过滤后图已保存至 results/QC_postfilter_violin.pdf\n")

# ============================================================================
# 6. 双细胞检测 (DoubletFinder)
# ============================================================================

cat("\n========== 双细胞检测 (DoubletFinder) ==========\n")
library(DoubletFinder)

run_doubletfinder <- function(seu, pK = NULL, pN = 0.25, doublet_rate = 0.075) {
  # 预处理
  seu <- NormalizeData(seu, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:10, verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)

  # 参数扫描
  sweep_res <- paramSweep(seu, PCs = 1:10, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  if (is.null(pK)) {
    pK <- as.numeric(as.character(sweep_stats$pK[which.max(sweep_stats$BCmetric)]))
  }
  cat("  Optimal pK:", pK, "\n")

  # 同型双细胞比例估计
  homotypic_prop <- modelHomotypic(seu$seurat_clusters)
  nExp_poi <- round(doublet_rate * ncol(seu))
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  cat("  Expected doublets:", nExp_poi, "| Adjusted:", nExp_poi_adj, "\n")

  # 运行 DoubletFinder
  seu <- doubletFinder_v3(seu,
    PCs = 1:10, pN = pN, pK = pK,
    nExp = nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE
  )

  # 提取分类列
  df_col <- grep("^DF\\.classifications", colnames(seu@meta.data), value = TRUE)[1]
  seu$doublet_class <- seu@meta.data[[df_col]]
  seu$doublet_class <- ifelse(seu$doublet_class == "Doublet", "Doublet", "Singlet")

  cat("  Singlets:", sum(seu$doublet_class == "Singlet"),
      "| Doublets:", sum(seu$doublet_class == "Doublet"), "\n")

  return(seu)
}

set.seed(42)
seu_cf   <- run_doubletfinder(seu_cf,   doublet_rate = 0.06)
seu_sham <- run_doubletfinder(seu_sham, doublet_rate = 0.06)
seu_acp  <- run_doubletfinder(seu_acp,  doublet_rate = 0.06)

# 去除双细胞
seu_cf   <- subset(seu_cf,   doublet_class == "Singlet")
seu_sham <- subset(seu_sham, doublet_class == "Singlet")
seu_acp  <- subset(seu_acp,  doublet_class == "Singlet")

cat("\n  DoubletFinder 完成，移除双细胞后:\n")
cat("  Mouse CF:", ncol(seu_cf), "cells\n")
cat("  Mouse Sham:", ncol(seu_sham), "cells\n")
cat("  Human ACP:", ncol(seu_acp), "cells\n")

# ============================================================================
# 7. 保存处理后的对象
# ============================================================================

cat("\n========== 保存数据 ==========\n")
dir.create(here("GSE206861/data/processed"), recursive = TRUE, showWarnings = FALSE)
saveRDS(seu_cf,   here("GSE206861/data/processed/seu_cf_qc.rds"))
saveRDS(seu_sham, here("GSE206861/data/processed/seu_sham_qc.rds"))
saveRDS(seu_acp,  here("GSE206861/data/processed/seu_acp_qc.rds"))

cat("\n✓ 数据读入与 QC 完成！\n")
cat("  保存: data/processed/seu_*_qc.rds\n")
cat("  接下来运行: 03_normalization_integration.R\n")
