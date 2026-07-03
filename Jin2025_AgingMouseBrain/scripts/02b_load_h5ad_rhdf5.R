# =============================================================================
# 用rhdf5直接读h5ad构建Seurat对象 (绕过SeuratDisk兼容性问题)
# =============================================================================
# h5ad 底层是 HDF5，rhdf5 可以直接读取
# AnnData结构: /X (稀疏矩阵), /obs (细胞元数据), /var (基因元数据)
# =============================================================================

load_h5ad_via_rhdf5 <- function(h5ad_path) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("rhdf5", update = FALSE, ask = FALSE)
  }
  library(rhdf5)
  library(Matrix)
  library(Seurat)
  library(data.table)

  cat(sprintf("  打开: %s\n", h5ad_path))

  # 读取细胞barcodes
  obs_keys <- h5read(h5ad_path, "obs/cell_label")
  cat(sprintf("  barcodes[1:3]: %s\n", paste(head(obs_keys, 3), collapse=", ")))

  # 读取基因名 (尝试多个可能路径)
  var_keys <- tryCatch(
    h5read(h5ad_path, "var/gene_identifier"),
    error = function(e) tryCatch(
      h5read(h5ad_path, "var/gene_symbol/categories"),
      error = function(e) paste0("gene_", seq_len(h5read(h5ad_path, "X/indptr") %>% length() - 1))
    )
  )

  cat(sprintf("  细胞: %d, 基因: %d\n", length(obs_keys), length(var_keys)))

  # 读取 X 稀疏矩阵 (CSR格式)
  # CSR: X/data (非零值), X/indices (列索引), X/indptr (行指针)
  X_data <- h5read(h5ad_path, "X/data")
  X_indices <- h5read(h5ad_path, "X/indices")
  X_indptr <- h5read(h5ad_path, "X/indptr")

  n_cells <- length(obs_keys)
  n_genes <- length(var_keys)

  cat(sprintf("  非零值: %d\n", length(X_data)))

  # 构建dgCMatrix (R的标准稀疏格式)
  # CSR: 第i行: data[indptr[i]:(indptr[i+1]-1)], cols=indices[indptr[i]:(indptr[i+1]-1)]
  # CSC: 第j列: 需要转换

  # 从CSR构建稀疏矩阵 — 使用Matrix::sparseMatrix
  # i = 行索引 (repeated for each non-zero)
  # j = 列索引 (from indices)
  # x = 值 (from data)
  # dims = c(n_cells, n_genes)

  cat("  构建稀疏矩阵...\n")

  # X_indptr: 长度 n_cells+1, CSR行指针 (每行的起止位置)
  # 检查indptr最后一个元素是否等于总非零数
  cat(sprintf("  X_indptr[1]=%d, X_indptr[last]=%d, length(X_data)=%d\n",
              X_indptr[1], X_indptr[length(X_indptr)], length(X_data)))

  # 确保indptr[last] == length(X_data)
  if (X_indptr[length(X_indptr)] != length(X_data)) {
    cat("  ⚠️ indptr[last] != length(data), 修正...\n")
    X_indptr <- c(X_indptr, length(X_data))
  }

  # 每行非零数 = diff(indptr)
  nnz_per_row <- diff(X_indptr)
  cat(sprintf("  nnz_per_row range: [%d, %d], sum=%d\n",
              min(nnz_per_row), max(nnz_per_row), sum(nnz_per_row)))

  # 检查负值
  if (any(nnz_per_row < 0)) {
    cat("  ❌ indptr不是非递减! 检查X_indptr...\n")
    print(head(which(nnz_per_row < 0), 10))
    stop("X_indptr must be non-decreasing")
  }

  row_vec <- rep(seq_len(n_cells), times = nnz_per_row)

  # 列索引: indices + 1 (Python 0-based -> R 1-based)
  col_vec <- as.integer(X_indices + 1)

  # 确保数据类型正确
  row_vec <- as.integer(row_vec)
  X_data_num <- as.numeric(X_data)

  cat(sprintf("  构建中... (i=%d, j=%d, x=%d elements)\n",
              length(row_vec), length(col_vec), length(X_data_num)))

  # 构建稀疏矩阵
  exp_mat <- sparseMatrix(
    i = row_vec,
    j = col_vec,
    x = X_data_num,
    dims = c(n_cells, n_genes),
    dimnames = list(obs_keys, var_keys)
  )

  # 清理临时变量
  rm(row_vec, col_vec, X_data_num, X_data, X_indices, X_indptr); gc()

  cat(sprintf("  矩阵: %d x %d, 非零率: %.2f%%\n",
              nrow(exp_mat), ncol(exp_mat),
              length(X_data) / (n_cells * n_genes) * 100))

  # 创建Seurat对象
  cat("  创建Seurat对象...\n")
  obj <- CreateSeuratObject(
    counts = exp_mat,
    assay = "RNA",
    meta.data = NULL
  )

  # 数据已经是log2标准化的，放入data slot
  obj[["RNA"]]$data <- obj[["RNA"]]$counts
  # 原始counts设为同样值(因为这是log2数据)
  obj[["RNA"]]$counts <- expm1(obj[["RNA"]]$counts)

  # 读取obs metadata (处理categorical编码)
  cat("  加载obs metadata...\n")
  obs_info <- h5ls(h5ad_path, recursive = FALSE)
  obs_datasets <- h5ls(paste0(h5ad_path, "/obs"), recursive = FALSE)

  for (ds in obs_datasets$name) {
    ds_path <- paste0("obs/", ds)
    info <- h5ls(h5ad_path, recursive = FALSE)

    # 检查是否为categorical (有categories和codes子项)
    cat_path <- paste0(h5ad_path, "/obs/", ds, "/categories")
    code_path <- paste0(h5ad_path, "/obs/", ds, "/codes")

    if (H5Lexists(h5ad_path, cat_path)) {
      # Categorical column: decode from codes
      cats <- h5read(h5ad_path, cat_path)
      codes <- h5read(h5ad_path, code_path)
      if (length(codes) == ncol(obj)) {
        obj[[ds]] <- cats[codes + 1]  # Python 0-based codes
      }
    } else {
      # Direct dataset
      vals <- tryCatch(
        h5read(h5ad_path, ds_path),
        error = function(e) NULL
      )
      if (!is.null(vals) && length(vals) == ncol(obj)) {
        obj[[ds]] <- vals
      }
    }
  }

  # 确保cell_type_major和age_group存在
  if (!"cell_type_major" %in% colnames(obj@meta.data) &&
      "cluster_name" %in% colnames(obj@meta.data)) {
    cat("  推断cell_type_major...\n")
    obj$cell_type_major <- sapply(obj$cluster_name, function(cn) {
      if (!is.character(cn)) return("Unclassified")
      if (grepl("Microglia", cn)) return("Microglia")
      if (grepl("Astro", cn) && !grepl("ependymal", cn, ignore.case=TRUE)) return("Astrocyte")
      if (grepl("Tanycyte", cn)) return("Tanycyte")
      if (grepl("Ependymal", cn)) return("Ependymal")
      if (grepl("NFOL", cn)) return("NFOL")
      if (grepl("MFOL", cn)) return("MFOL")
      if (grepl("MOL", cn)) return("MOL")
      if (grepl("OPC", cn)) return("OPC")
      if (grepl("COP", cn)) return("COP")
      if (grepl("Endo", cn)) return("Endothelial")
      if (grepl("SMC", cn)) return("SMC")
      if (grepl("VLMC", cn)) return("VLMC")
      if (grepl("Peri", cn)) return("Pericyte")
      if (grepl("BAM", cn)) return("BAM")
      if (grepl("_DC|DC_", cn)) return("Dendritic_Cell")
      if (grepl("ABC", cn)) return("ABC")
      if (grepl("T cells", cn)) return("T_cell")
      if (grepl("Glut", cn) && !grepl("Gaba|GABA", cn)) return("Glutamatergic_Neuron")
      if (grepl("Gaba|GABA", cn)) return("GABAergic_Neuron")
      return("Other")
    })
  }

  cat(sprintf("  ✓ Seurat对象: %d cells, %d genes\n", ncol(obj), nrow(obj)))

  # 清理
  rm(exp_mat); gc()

  return(obj)
}

# =============================================================================
# 测试: 加载下丘脑数据
# =============================================================================
cat("============================================\n")
cat(" rhdf5直接读取h5ad\n")
cat("============================================\n\n")

h5ad_anat <- "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_anatomical_52696.h5ad"

hypo_anat <- load_h5ad_via_rhdf5(h5ad_anat)

# 保存为RDS
saveRDS(hypo_anat, "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_anatomical_seurat.rds")
cat(sprintf("\n✓ 已保存: hypo_anatomical_seurat.rds\n"))
cat(sprintf("  细胞: %d\n", ncol(hypo_anat)))
cat(sprintf("  基因: %d\n", nrow(hypo_anat)))
