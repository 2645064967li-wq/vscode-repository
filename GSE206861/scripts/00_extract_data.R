###############################################################################
# GSE206861 数据解压与整理
# 将下载好的 GSE206861_RAW.tar 解压并按样本分目录存放
#
# 使用方式:
#   在 R 中运行: source("scripts/00_extract_data.R")
#   或在终端运行: tar -xf GSE206861_RAW.tar -C data/
###############################################################################

# 检查是否有 tar 文件
tar_file <- "GSE206861_RAW.tar"
if (!file.exists(tar_file)) {
  stop("找不到 ", tar_file,
       "\n请先从浏览器下载: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206861/suppl/GSE206861_RAW.tar",
       "\n或使用下载工具下载后放在 GSE206861/ 目录下")
}

cat("找到:", tar_file, "(", round(file.info(tar_file)$size / 1024^2, 1), "MB )\n")

# 创建数据目录
dir.create("data", recursive = TRUE, showWarnings = FALSE)

# 解压到 data/ 目录
cat("正在解压...\n")
utils::untar(tar_file, exdir = "data")
cat("解压完成!\n")

# 列出解压后的文件
files <- list.files("data", pattern = "\\.gz$", recursive = TRUE)
cat("\n解压文件列表:\n")
for (f in files) {
  cat("  ", f, "(", round(file.info(file.path("data", f))$size / 1024^2, 1), "MB )\n")
}

# ============================================================================
# 按样本整理目录
# ============================================================================

cat("\n正在按样本整理目录...\n")

# 定义样本信息
samples <- list(
  "GSM6265811_Mouse_Cystic-fluid" = c(
    "GSM6265811_Mouse_Cystic-fluid_barcodes.tsv.gz",
    "GSM6265811_Mouse_Cystic-fluid_features.tsv.gz",
    "GSM6265811_Mouse_Cystic-fluid_matrix.mtx.gz"
  ),
  "GSM6265812_Mouse_Sham" = c(
    "GSM6265812_Mouse_Sham_barcodes.tsv.gz",
    "GSM6265812_Mouse_Sham_features.tsv.gz",
    "GSM6265812_Mouse_Sham_matrix.mtx.gz"
  ),
  "GSM6265813_Human_ACP" = c(
    "GSM6265813_Human_ACP_barcodes.tsv.gz",
    "GSM6265813_Human_ACP_features.tsv.gz",
    "GSM6265813_Human_ACP_matrix.mtx.gz"
  )
)

# 移动文件到子目录
for (sample_dir in names(samples)) {
  dir.create(file.path("data", sample_dir), recursive = TRUE, showWarnings = FALSE)

  for (f in samples[[sample_dir]]) {
    old_path <- file.path("data", f)
    new_path <- file.path("data", sample_dir, f)

    if (file.exists(old_path)) {
      # 去掉样本名前缀，重命名为标准 10x 格式
      simple_name <- sub(paste0(sample_dir, "_"), "", f)
      simple_path <- file.path("data", sample_dir, simple_name)

      file.rename(old_path, new_path)
      file.rename(new_path, simple_path)

      cat("  ", sample_dir, "/", simple_name, "\n")
    }
  }
}

# ============================================================================
# 检查 10x 文件完整性
# ============================================================================

cat("\n========== 文件完整性检查 ==========\n")

for (sample_dir in names(samples)) {
  cat("\n", sample_dir, ":\n")

  barcodes <- file.path("data", sample_dir, "barcodes.tsv.gz")
  features <- file.path("data", sample_dir, "features.tsv.gz")
  matrix   <- file.path("data", sample_dir, "matrix.mtx.gz")

  check_file <- function(path, label) {
    if (file.exists(path)) {
      cat("  [OK] ", label, "(", round(file.info(path)$size / 1024, 1), "KB )\n")
    } else {
      cat("  [MISSING] ", label, "\n")
    }
  }

  check_file(barcodes, "barcodes.tsv.gz")
  check_file(features, "features.tsv.gz")
  check_file(matrix, "matrix.mtx.gz")

  # 检查 matrix 维度
  if (file.exists(matrix)) {
    # 读取 matrix.mtx 前几行获取维度信息
    con <- gzfile(matrix, "r")
    header_lines <- readLines(con, n = 10)
    close(con)
    # 第三行是维度信息
    dim_line <- header_lines[3]
    cat("  Matrix dims:", dim_line, "\n")
  }
}

cat("\n✓ 数据整理完成！可以运行 01_install_packages.R 开始分析。\n")
