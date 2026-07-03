###############################################################################
# GSE206861 单细胞分析 — 环境搭建脚本
# 文献: ACP cyst fluid triggers microglia activation (PMID: 35525962)
# 分析: R 4.x + Seurat 5
###############################################################################

# ============================================================================
# Step 0: 设置镜像（国内用户建议用镜像加速）
# ============================================================================
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(timeout = 600)

# ============================================================================
# Step 1: 检查 R 版本
# ============================================================================
cat("\n========== R 版本信息 ==========\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")

required_r <- "4.1.0"
if (compareVersion(paste0(R.version$major, ".", R.version$minor), required_r) < 0) {
  stop("需要 R >= ", required_r, "，当前版本: ", R.version.string)
}

# ============================================================================
# Step 2: 安装 BiocManager（Bioconductor 包管理器）
# ============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("\n>>> 安装 BiocManager...\n")
  install.packages("BiocManager")
}

# ============================================================================
# Step 3: 安装 CRAN 包
# ============================================================================
cran_packages <- c(
  # 核心框架
  "Seurat",           # 单细胞分析核心
  "SeuratObject",     # Seurat 数据对象
  "SeuratWrappers",   # Seurat 扩展（Monocle3 转换等）

  # 数据处理与绘图
  "tidyverse",        # dplyr + ggplot2 + tidyr + readr + ...
  "patchwork",        # 图组合
  "cowplot",          # 出版级拼图
  "viridis",          # 色盲友好配色
  "ggrepel",          # 标签防重叠
  "ggrastr",          # 光栅化绘图（大数据集更流畅）
  "RColorBrewer",     # 配色

  # 分析工具
  "harmony",          # 批次校正 / 整合
  "future",           # 并行计算
  "Matrix",           # 稀疏矩阵
  "irlba",            # 快速 SVD

  # 辅助
  "here",             # 相对路径管理
  "qs",               # 快速序列化（保存/读取大对象）
  "devtools"          # 从 GitHub 安装包
)

cat("\n>>> 安装 CRAN 包（共", length(cran_packages), "个）...\n")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  安装:", pkg, "\n")
    tryCatch(
      install.packages(pkg),
      error = function(e) cat("  !! 安装失败:", pkg, "-", e$message, "\n")
    )
  } else {
    cat("  [已安装]", pkg, "\n")
  }
}

# ============================================================================
# Step 4: 安装 Bioconductor 包
# ============================================================================
bioc_packages <- c(
  # 自动细胞注释
  "SingleR",          # 基于参考数据集的自动注释
  "celldex",          # 参考表达数据集

  # 富集分析
  "clusterProfiler",  # GO/KEGG 富集
  "enrichplot",       # 富集结果可视化
  "DOSE",             # 疾病本体富集

  # 基因 ID 转换
  "org.Mm.eg.db",     # 小鼠基因
  "org.Hs.eg.db",     # 人类基因

  # 基因集
  "msigdbr",          # MSigDB 基因集（含 Hallmark, GO, KEGG 等）

  # 拟时序
  "slingshot",        # 轨迹推断

  # 其他
  "scater",           # 单细胞数据基础操作
  "scran",            # 单细胞归一化
  "DropletUtils",     # 液滴数据工具
  "GEOquery",         # GEO 数据下载
  "BiocParallel",     # 并行计算
  "limma"             # 差异分析（备用）
)

cat("\n>>> 安装 Bioconductor 包（共", length(bioc_packages), "个）...\n")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  安装:", pkg, "\n")
    tryCatch(
      BiocManager::install(pkg, update = FALSE, ask = FALSE),
      error = function(e) cat("  !! 安装失败:", pkg, "-", e$message, "\n")
    )
  } else {
    cat("  [已安装]", pkg, "\n")
  }
}

# ============================================================================
# Step 5: 从 GitHub 安装（如果需要）
# ============================================================================
github_packages <- list(
  # DoubletFinder — 双细胞检测（CRAN 没有）
  "DoubletFinder" = "chris-mcginnis-uoft/DoubletFinder",
  # CellChat — 细胞通讯分析
  "CellChat"       = "jinworks/CellChat"
)

cat("\n>>> 安装 GitHub 包（共", length(github_packages), "个）...\n")

for (pkg_name in names(github_packages)) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat("  安装:", pkg_name, "\n")
    tryCatch(
      devtools::install_github(github_packages[[pkg_name]], upgrade = "never"),
      error = function(e) cat("  !! 安装失败:", pkg_name, "-", e$message, "\n")
    )
  } else {
    cat("  [已安装]", pkg_name, "\n")
  }
}

# ============================================================================
# Step 6: 安装 Monocle3（从 Bioconductor，注意特定依赖）
# ============================================================================
cat("\n>>> 安装 Monocle3...\n")
if (!requireNamespace("monocle3", quietly = TRUE)) {
  tryCatch(
    {
      BiocManager::install(c("DelayedArray", "DelayedMatrixStats",
                              "SingleCellExperiment", "SummarizedExperiment",
                              "batchelor", "HDF5Array", "terra",
                              "ggrastr", "BiocGenerics"),
                            update = FALSE, ask = FALSE)
      devtools::install_github("cole-trapnell-lab/monocle3", upgrade = "never")
    },
    error = function(e) cat("  !! Monocle3 安装失败:", e$message, "\n")
  )
} else {
  cat("  [已安装] monocle3\n")
}

# ============================================================================
# Step 7: 验证所有包
# ============================================================================
cat("\n========== 包验证 ==========\n")

all_check_packages <- c(
  cran_packages,
  bioc_packages,
  names(github_packages),
  "monocle3"
)

# 去除不存在的包名
all_check_packages <- unique(all_check_packages)

failures <- c()
for (pkg in all_check_packages) {
  loaded <- suppressPackageStartupMessages(
    requireNamespace(pkg, quietly = TRUE)
  )
  if (loaded) {
    cat("  [OK]", pkg, "\n")
  } else {
    cat("  [FAIL]", pkg, "\n")
    failures <- c(failures, pkg)
  }
}

if (length(failures) > 0) {
  cat("\n⚠ 以下包未成功安装:", paste(failures, collapse = ", "), "\n")
  cat("请手动检查安装错误信息，可能需要额外的系统依赖。\n")
  cat("参考: https://satijalab.org/seurat/articles/install.html\n")
} else {
  cat("\n✓ 所有包安装成功！可以开始分析了。\n")
}

# ============================================================================
# Step 8: 设置并行计算
# ============================================================================
cat("\n========== 并行设置 ==========\n")
cpu_cores <- parallel::detectCores()
cat("检测到 CPU 核心数:", cpu_cores, "\n")
cat("建议使用:", max(1, cpu_cores - 2), "个核心用于分析\n")

# ============================================================================
# 完成
# ============================================================================
cat("\n========== 环境搭建完成 ==========\n")
cat("如果所有包都显示 [OK]，请运行下一个脚本: 02_load_data_QC.R\n")
