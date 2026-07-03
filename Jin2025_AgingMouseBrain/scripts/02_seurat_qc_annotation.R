# =============================================================================
# Phase 2: 下丘脑Seurat QC与注释
# =============================================================================
# 处理两个label集合: anatomical(52K) 和 dissection(26K)
# =============================================================================
set.seed(42)

# ---- 0. 环境 ----
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_cran <- c("Seurat", "dplyr", "ggplot2", "patchwork", "data.table",
               "Matrix", "cowplot", "scales")
pkgs_bioc <- c("SingleCellExperiment")
pkgs_gh <- c("mojaveazure/seurat-disk")

for (pkg in pkgs_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
for (pkg in pkgs_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update=FALSE, ask=FALSE)
}
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("mojaveazure/seurat-disk", upgrade="never")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
})

cat("\n============================================\n")
cat(" Phase 2: 下丘脑 Seurat QC + 注释\n")
cat("============================================\n\n")

# ---- 1. 路径 ----
DATA_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus"
RESULT_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/results"
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 2. 处理函数 ----
process_label_set <- function(label_name, h5ad_path) {
  cat(sprintf("\n═══════════════════════════════════════\n"))
  cat(sprintf(" 处理: %s\n", label_name))
  cat(sprintf("═══════════════════════════════════════\n\n"))

  # 2a. 加载h5ad
  cat("[Step 1] 加载h5ad...\n")
  stopifnot(file.exists(h5ad_path))
  h5seurat_file <- sub("\\.h5ad$", ".h5seurat", h5ad_path)

  if (!file.exists(h5seurat_file)) {
    cat("  转换 h5ad -> h5Seurat...\n")
    Convert(h5ad_path, dest = h5seurat_file, overwrite = TRUE)
  }
  obj <- LoadH5Seurat(h5seurat_file)
  cat(sprintf("  细胞数: %s\n", format(ncol(obj), big.mark=",")))
  cat(sprintf("  基因数: %s\n", format(nrow(obj), big.mark=",")))

  # 2b. 检查meta.data中的注释
  cat("\n[Step 2] 注释检查...\n")

  # 确保有关键列
  required_cols <- c("cell_type_major", "age_group", "sex_group", "group")
  for (col in required_cols) {
    if (col %in% colnames(obj@meta.data)) {
      cat(sprintf("  ✓ %s: %s unique values\n", col,
                  paste(unique(obj@meta.data[[col]]), collapse=", ")))
    } else {
      cat(sprintf("  ⚠️ %s: MISSING\n", col))
    }
  }

  # 过滤掉age_group=Unknown的细胞
  n_before <- ncol(obj)
  if ("age_group" %in% colnames(obj@meta.data)) {
    obj <- subset(obj, age_group != "Unknown")
    cat(sprintf("  过滤Unknown age: %d -> %d cells\n", n_before, ncol(obj)))
  }

  # 2c. QC指标
  cat("\n[Step 3] QC指标...\n")

  # 计算线粒体比例
  mito_genes <- grep("^mt-|^Mt-", rownames(obj), value = TRUE)
  if (length(mito_genes) > 0) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mito_genes[1:min(100, length(mito_genes))])
  } else {
    obj[["percent.mt"]] <- 0
    cat("  ⚠️ 未找到线粒体基因\n")
  }

  # QC分布 (不硬过滤,因为已经是高质量数据)
  cat("  QC指标分布:\n")
  for (metric in c("gene_count", "umi_count", "doublet_score")) {
    if (metric %in% colnames(obj@meta.data)) {
      vals <- obj@meta.data[[metric]]
      cat(sprintf("    %s: median=%.1f, IQR=[%.1f-%.1f]\n",
                  metric, median(vals, na.rm=TRUE),
                  quantile(vals, 0.25, na.rm=TRUE),
                  quantile(vals, 0.75, na.rm=TRUE)))
    }
  }

  # 2d. 标准化和降维
  cat("\n[Step 4] 标准化 + 降维...\n")

  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000)

  # 回归性别和UMI
  vars_to_regress <- "umi_count"
  if ("umi_count" %in% colnames(obj@meta.data)) {
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
    cat(sprintf("  高变基因: %d\n", length(VariableFeatures(obj))))
    obj <- ScaleData(obj, vars.to.regress = vars_to_regress)
  } else {
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
    cat(sprintf("  高变基因: %d\n", length(VariableFeatures(obj))))
    obj <- ScaleData(obj)
  }

  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50)
  obj <- RunUMAP(obj, dims = 1:30, n.neighbors = 30, min.dist = 0.3)

  # 2e. 按cell_type_major统计
  cat("\n[Step 5] 细胞类型汇总...\n")

  # 确保cell_type_major存在
  if (!"cell_type_major" %in% colnames(obj@meta.data)) {
    # 从cluster_name推断
    if ("cluster_name" %in% colnames(obj@meta.data)) {
      cat("  从cluster_name推断cell_type_major...\n")
      obj$cell_type_major <- sapply(obj$cluster_name, function(cn) {
        if (!is.character(cn)) return("Unclassified")
        if (grepl("Microglia", cn)) return("Microglia")
        if (grepl("Astro", cn) && !grepl("ependymal", cn, ignore.case=TRUE))
          return("Astrocyte")
        if (grepl("Tanycyte", cn)) return("Tanycyte")
        if (grepl("Ependymal", cn)) return("Ependymal")
        if (grepl("NFOL", cn)) return("NFOL")
        if (grepl("MFOL", cn)) return("MFOL")
        if (grepl("MOL", cn)) return("MOL")
        if (grepl("OPC", cn)) return("OPC")
        if (grepl("COP", cn)) return("COP")
        if (grepl("Endo", cn)) return("Endothelial")
        if (grepl("SMC|SMC-", cn)) return("SMC")
        if (grepl("VLMC", cn)) return("VLMC")
        if (grepl("Peri", cn)) return("Pericyte")
        if (grepl("BAM", cn)) return("BAM")
        if (grepl("_DC|DC_", cn)) return("Dendritic_Cell")
        if (grepl("ABC", cn)) return("ABC")
        if (grepl("T cells", cn)) return("T_cell")
        if (grepl("Glut", cn) && !grepl("Gaba|GABA", cn))
          return("Glutamatergic_Neuron")
        if (grepl("Gaba|GABA", cn)) return("GABAergic_Neuron")
        return("Other")
      })
    }
  }

  # 生成汇总表
  Idents(obj) <- "cell_type_major"
  ct_summary <- obj@meta.data %>%
    group_by(cell_type_major) %>%
    summarise(
      total_cells = n(),
      adult_cells = sum(age_group == "Adult(2m)", na.rm = TRUE),
      aged_cells  = sum(age_group == "Aged(18m)", na.rm = TRUE),
      female_cells = sum(sex_group == "Female", na.rm = TRUE),
      male_cells   = sum(sex_group == "Male", na.rm = TRUE),
      aged_pct = round(aged_cells / total_cells * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(total_cells))

  print(ct_summary, n = 30)

  # 保存汇总
  fwrite(ct_summary, file.path(RESULT_DIR,
        sprintf("cell_type_summary_%s.csv", label_name)))

  # 2f. 快速UMAP可视化
  cat("\n[Step 6] 生成基础UMAP...\n")

  p1 <- DimPlot(obj, group.by = "cell_type_major", label = TRUE,
                repel = TRUE, pt.size = 0.3) +
    ggtitle(sprintf("Hypothalamus: Cell Types (%s)", label_name)) +
    theme_minimal(base_size = 10) +
    theme(legend.text = element_text(size = 7))

  p2 <- DimPlot(obj, group.by = "age_group",
                cols = c("Adult(2m)" = "#4DBBD5", "Aged(18m)" = "#E64B35"),
                pt.size = 0.3) +
    ggtitle(sprintf("Hypothalamus: Age (%s)", label_name)) +
    theme_minimal(base_size = 10)

  ggsave(file.path(RESULT_DIR, sprintf("umap_celltype_%s.png", label_name)),
         p1, width = 12, height = 8, dpi = 150)
  ggsave(file.path(RESULT_DIR, sprintf("umap_age_%s.png", label_name)),
         p2, width = 8, height = 7, dpi = 150)
  cat("  ✓ UMAP已保存\n")

  # 2g. 保存Seurat对象
  cat("\n[Step 7] 保存Seurat RDS...\n")
  rds_file <- file.path(DATA_DIR, sprintf("hypo_%s_seurat.rds", label_name))
  saveRDS(obj, rds_file)
  cat(sprintf("  ✓ 已保存: %s\n", rds_file))

  # 清理
  gc()
  return(obj)
}

# ---- 3. 运行 ----
# 处理 anatomical (52K)
h5ad_anat <- file.path(DATA_DIR, "hypo_anatomical_52696.h5ad")
if (file.exists(h5ad_anat)) {
  hypo_anat <- process_label_set("anatomical", h5ad_anat)
}

# 处理 dissection (26K)
h5ad_diss <- file.path(DATA_DIR, "hypo_dissection_26105.h5ad")
if (file.exists(h5ad_diss)) {
  hypo_diss <- process_label_set("dissection", h5ad_diss)
}

cat("\n============================================\n")
cat(" Phase 2 完成!\n")
cat("============================================\n")
