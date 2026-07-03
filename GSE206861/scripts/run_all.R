###############################################################################
# GSE206861 单细胞分析 — 主运行脚本
# 一键执行完整分析流程
#
# 使用方法:
#   1. 先在 RStudio 中运行 01_install_packages.R 安装依赖
#   2. 下载 GSE206861_RAW.tar 并运行 00_extract_data.R 解压
#   3. 运行本脚本: source("scripts/run_all.R")
#     或者分步运行各个脚本
#
# 文献: ACP cyst fluid triggers microglia activation (PMID: 35525962)
# 数据集: GEO GSE206861
###############################################################################

library(here)

# ============================================================================
# 配置参数（可根据需要修改）
# ============================================================================

# QC 参数
QC_PARAMS <- list(
  mouse_nFeature_low  = 200,
  mouse_nFeature_high = 6000,
  mouse_nCount_low    = 500,
  mouse_nCount_high   = 50000,
  mouse_percent_mt    = 20,

  human_nFeature_low  = 200,
  human_nFeature_high = 5000,
  human_nCount_low    = 500,
  human_nCount_high   = 40000,
  human_percent_mt    = 20
)

# 聚类参数
CLUSTER_PARAMS <- list(
  mouse_dims       = 30,    # Harmony/PCA 维度
  mouse_resolution = 0.6,   # 聚类分辨率
  human_dims       = 20,
  human_resolution = 0.6
)

# 差异分析参数
DE_PARAMS <- list(
  min_pct       = 0.1,
  logfc_threshold = 0.25,
  test_method   = "wilcox"
)

# ============================================================================
# 设置
# ============================================================================
setwd(here("GSE206861"))
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

cat("=============================================\n")
cat("  GSE206861 单细胞分析流程\n")
cat("  PMID: 35525962\n")
cat("=============================================\n\n")

cat("QC 参数:\n")
print(QC_PARAMS)
cat("\n聚类参数:\n")
print(CLUSTER_PARAMS)

# ============================================================================
# Step 0: 数据解压（如果尚未解压）
# ============================================================================
cat("\n\n========== Step 0: 检查数据 ==========\n")

data_ready <- TRUE
required_files <- c(
  "data/GSM6265811_Mouse_Cystic-fluid/barcodes.tsv.gz",
  "data/GSM6265811_Mouse_Cystic-fluid/features.tsv.gz",
  "data/GSM6265811_Mouse_Cystic-fluid/matrix.mtx.gz",
  "data/GSM6265812_Mouse_Sham/barcodes.tsv.gz",
  "data/GSM6265812_Mouse_Sham/features.tsv.gz",
  "data/GSM6265812_Mouse_Sham/matrix.mtx.gz",
  "data/GSM6265813_Human_ACP/barcodes.tsv.gz",
  "data/GSM6265813_Human_ACP/features.tsv.gz",
  "data/GSM6265813_Human_ACP/matrix.mtx.gz"
)

for (f in required_files) {
  if (!file.exists(f)) {
    cat("  [MISSING]", f, "\n")
    data_ready <- FALSE
  }
}

if (!data_ready) {
  cat("\n数据文件不完整。请先运行: source('scripts/00_extract_data.R')\n")
  cat("如果尚未下载 tar 文件，请从浏览器下载:\n")
  cat("  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206861/suppl/GSE206861_RAW.tar\n")
  stop("数据准备不完整。")
} else {
  cat("  数据文件完整，继续分析。\n")
}

# ============================================================================
# Step 1: 数据读入 + QC
# ============================================================================
cat("\n\n========== Step 1: 数据读入 + QC ==========\n")
source("scripts/02_load_data_QC.R")

# ============================================================================
# Step 2: 标准化 + 整合 + 降维 + 聚类
# ============================================================================
cat("\n\n========== Step 2: 标准化 + 整合 + 聚类 ==========\n")
source("scripts/03_normalization_integration.R")

# ============================================================================
# Step 3: 细胞注释
# ============================================================================
cat("\n\n========== Step 3: 细胞注释 ==========\n")
source("scripts/04_cell_annotation.R")

# ============================================================================
# Step 4: 差异分析 + 富集分析
# ============================================================================
cat("\n\n========== Step 4: 差异分析 + 富集 ==========\n")
source("scripts/05_differential_expression_enrichment.R")

# ============================================================================
# Step 5: 进阶分析（可选，耗时较长）
# ============================================================================
cat("\n\n========== Step 5: 进阶分析 ==========\n")

run_advanced <- TRUE  # 设置为 FALSE 跳过拟时序和细胞通讯分析

if (run_advanced) {
  cat("正在运行拟时序和细胞通讯分析...\n")
  source("scripts/06_advanced_analysis.R")
} else {
  cat("跳过进阶分析。设置 run_advanced = TRUE 以运行。\n")
}

# ============================================================================
# 完成
# ============================================================================
cat("\n=============================================\n")
cat("  全部分析完成！\n")
cat("  结果文件位于: results/\n")
cat("  处理数据位于: data/processed/\n")
cat("=============================================\n")

# 列出所有结果文件
cat("\n生成的结果文件:\n")
result_files <- list.files("results", recursive = TRUE, full.names = TRUE)
for (f in sort(result_files)) {
  cat("  ", f, "\n")
}
