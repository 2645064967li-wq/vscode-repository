# =============================================================================
# R 版本: 下载 Jin et al. 2025 衰老小鼠全脑数据 - Metadata
# =============================================================================
# 使用方法:
#   Rscript scripts/01_download_metadata.R
#   或在 RStudio 中: source("scripts/01_download_metadata.R")
# =============================================================================

# 配置 ----
BASE_URL <- "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
PROJECT_DIR <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
if (PROJECT_DIR == "" || is.null(PROJECT_DIR)) {
  PROJECT_DIR <- getwd()
}
DATA_DIR <- file.path(PROJECT_DIR, "data")

# ============================================================================
# 辅助函数 ----
# ============================================================================
download_with_retry <- function(url, destfile, max_retries = 3) {
  for (i in 1:max_retries) {
    tryCatch({
      download.file(url, destfile, mode = "wb", quiet = FALSE)
      return(TRUE)
    }, error = function(e) {
      cat(sprintf("  重试 %d/%d: %s\n", i, max_retries, e$message))
      Sys.sleep(5)
    })
  }
  cat(sprintf("  ✗ 下载失败 (重试%d次后): %s\n", max_retries, basename(destfile)))
  return(FALSE)
}

# ============================================================================
# 1. 单细胞转录组 Metadata (Zeng-Aging-Mouse-10Xv3) ----
# ============================================================================
cat("\n[1/3] 下载单细胞 metadata (Zeng-Aging-Mouse-10Xv3)...\n")
cat("      版本: 20250131\n")

sc_meta_dir <- file.path(DATA_DIR, "single_cell", "metadata")
dir.create(sc_meta_dir, recursive = TRUE, showWarnings = FALSE)

sc_files <- c(
  "cell_metadata.csv",
  "cell_annotation_colors.csv",
  "cell_cluster_annotations.csv",
  "cluster.csv",
  "donor.csv",
  "library.csv",
  "value_sets.csv",
  "views/example_genes_all_cells_expression.csv"
)

sc_meta_url <- file.path(BASE_URL, "metadata", "Zeng-Aging-Mouse-10Xv3", "20250131")

for (file in sc_files) {
  cat(sprintf("  -> 下载 %s...\n", file))
  dest <- file.path(sc_meta_dir, file)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  download_with_retry(file.path(sc_meta_url, file), dest)
}

cat("  ✓ 单细胞 metadata 完成!\n")

# ============================================================================
# 2. 单细胞 Taxonomy Metadata ----
# ============================================================================
cat("\n[2/3] 下载 Taxonomy metadata (WMB-taxonomy)...\n")
cat("      版本: 20241130\n")

tax_meta_dir <- file.path(DATA_DIR, "single_cell", "taxonomy")
dir.create(tax_meta_dir, recursive = TRUE, showWarnings = FALSE)

tax_files <- c(
  "aging_degenes.csv",
  "cell_cluster_mapping_annotations.csv",
  "cell_cross_mapping_annotations.csv",
  "cluster_mapping.csv",
  "cluster_mapping_pivot.csv"
)

tax_url <- file.path(BASE_URL, "metadata", "Zeng-Aging-Mouse-WMB-taxonomy", "20241130")

for (file in tax_files) {
  cat(sprintf("  -> 下载 %s...\n", file))
  dest <- file.path(tax_meta_dir, file)
  download_with_retry(file.path(tax_url, file), dest)
}

cat("  ✓ Taxonomy metadata 完成!\n")

# ============================================================================
# 3. MERFISH 空间转录组 Metadata ----
# ============================================================================
cat("\n[3/3] 下载 MERFISH 空间转录组 metadata...\n")

merfish_meta_dir <- file.path(DATA_DIR, "spatial", "MERFISH_638850", "metadata")
dir.create(merfish_meta_dir, recursive = TRUE, showWarnings = FALSE)

# Cell metadata
merfish_url <- file.path(BASE_URL, "metadata", "MERFISH-C57BL6J-638850", "20241115")
merfish_files <- c("cell_metadata.csv", "gene.csv",
                   "views/cell_metadata_with_cluster_annotation.csv")

for (file in merfish_files) {
  cat(sprintf("  -> 下载 %s...\n", file))
  dest <- file.path(merfish_meta_dir, file)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  download_with_retry(file.path(merfish_url, file), dest)
}

# CCF coordinates
cat("  -> 下载 ccf_coordinates.csv (CCF空间坐标)...\n")
ccf_url <- file.path(BASE_URL, "metadata", "MERFISH-C57BL6J-638850-CCF", "20231215")
download_with_retry(
  file.path(ccf_url, "ccf_coordinates.csv"),
  file.path(merfish_meta_dir, "ccf_coordinates.csv")
)

# Imputed gene list
cat("  -> 下载 imputed gene list...\n")
imp_url <- file.path(BASE_URL, "metadata", "MERFISH-C57BL6J-638850-imputed", "20240831")
download_with_retry(
  file.path(imp_url, "gene.csv"),
  file.path(merfish_meta_dir, "gene_imputed.csv")
)

cat("  ✓ MERFISH metadata 完成!\n")

# ============================================================================
# 4. 快速检查 ----
# ============================================================================
cat("\n============================================\n")
cat(" 下载完成! 快速数据检查:\n")
cat("============================================\n\n")

# 检查单细胞 metadata
cell_meta <- file.path(sc_meta_dir, "cell_metadata.csv")
if (file.exists(cell_meta)) {
  cat(sprintf("cell_metadata.csv: %s\n", file.info(cell_meta)$size))
  # 快速看一眼前几列
  tryCatch({
    meta_head <- read.csv(cell_meta, nrows = 5)
    cat(sprintf("  列数: %d\n", ncol(meta_head)))
    cat(sprintf("  部分列名: %s...\n",
                paste(head(colnames(meta_head), 10), collapse = ", ")))
  }, error = function(e) cat(sprintf("  读取错误: %s\n", e$message)))
} else {
  cat("  ✗ cell_metadata.csv 未找到!\n")
}

# 检查 aging_degenes
deg_file <- file.path(tax_meta_dir, "aging_degenes.csv")
if (file.exists(deg_file)) {
  cat(sprintf("\naging_degenes.csv: %s\n", file.info(deg_file)$size))
  tryCatch({
    deg <- read.csv(deg_file, nrows = 5)
    cat(sprintf("  年龄差异基因总数: 预计2,449个\n"))
  }, error = function(e) invisible())
}

# 检查 MERFISH
ccf_file <- file.path(merfish_meta_dir, "ccf_coordinates.csv")
if (file.exists(ccf_file)) {
  cat(sprintf("\nccf_coordinates.csv: %s (%s)\n",
              file.info(ccf_file)$size, "CCFv3空间坐标"))
}

cat("\n")
cat("下一步:\n")
cat("  1. 运行 scripts/02_download_expression_matrices.sh 下载表达矩阵\n")
cat("  2. 或运行 Rscript scripts/filter_hypothalamus.R 提取下丘脑数据\n")
