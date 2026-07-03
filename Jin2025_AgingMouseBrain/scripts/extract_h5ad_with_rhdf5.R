# =============================================================================
# 用rhdf5从h5ad中提取下丘脑Microglia细胞 (纯R, 无需Python)
# =============================================================================
# h5ad底层是HDF5格式, rhdf5可以直接读
# 只提取2,047个cells, 避免整数溢出
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# 安装rhdf5
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("rhdf5", quietly = TRUE))
  BiocManager::install("rhdf5", update = FALSE, ask = FALSE)

library(rhdf5)

H5AD_FILE <- "d:/decrepitude mouse hypothamulas/Zeng-Aging-Mouse-10Xv3-log2.h5ad"
BARCODE_FILE <- "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_microglia_barcodes.txt"
OUT_RDS <- "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_microglia_raw_counts.rds"
OUT_META <- "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_microglia_meta.rds"

cat("============================================\n")
cat(" rhdf5 提取下丘脑Microglia表达数据\n")
cat("============================================\n\n")

# 1. 读入目标barcodes
target_bc <- readLines(BARCODE_FILE)
cat(sprintf("目标cells: %d\n", length(target_bc)))

# 2. 打开h5ad, 读取cell barcodes
cat("\n[1] 读取h5ad cell index...\n")
h5ls(H5AD_FILE, recursive = FALSE)  # 显示顶层结构

# obs是cell metadata
all_barcodes <- h5read(H5AD_FILE, "obs/_index")
cat(sprintf("总cells: %d\n", length(all_barcodes)))

# 或者从obs的index列读取
if (length(all_barcodes) == 0) {
  # 尝试其他位置
  obs_names <- h5read(H5AD_FILE, "obs/cell_label")
  if (length(obs_names) > 0) all_barcodes <- obs_names
}
cat(sprintf("  barcodes[1:5]: %s\n", paste(head(all_barcodes, 5), collapse = ", ")))
cat(sprintf("  target[1:3]: %s\n", paste(head(target_bc, 3), collapse = ", ")))

# 3. 匹配barcodes, 找到目标细胞的索引
cat("\n[2] 匹配barcodes...\n")
matched_idx <- which(all_barcodes %in% target_bc)
cat(sprintf("  匹配到: %d cells\n", length(matched_idx)))

if (length(matched_idx) == 0) {
  cat("\n  ⚠️ 精确匹配失败, 尝试模糊匹配...\n")
  # 去掉可能的library后缀
  target_stems <- gsub("-.*$", "", target_bc)
  all_stems <- gsub("-.*$", "", all_barcodes)
  matched_idx <- which(all_stems %in% target_stems)
  cat(sprintf("  模糊匹配: %d cells\n", length(matched_idx)))
}

if (length(matched_idx) == 0) {
  cat("\n  ❌ 无法匹配! 可能需要不同的barcode格式\n")
  quit(status = 1)
}

# 4. 读取稀疏矩阵 X (CSR格式)
cat("\n[3] 读取X稀疏矩阵 (只提取目标cells)...\n")

# AnnData的X在CSR格式: /X/data, /X/indices, /X/indptr
X_data <- h5read(H5AD_FILE, "X/data")
X_indices <- h5read(H5AD_FILE, "X/indices")
X_indptr <- h5read(H5AD_FILE, "X/indptr")

cat(sprintf("  X/data length: %d\n", length(X_data)))
cat(sprintf("  X/indices length: %d\n", length(X_indices)))
cat(sprintf("  X/indptr length: %d (should be n_cells+1=%d)\n",
            length(X_indptr), length(all_barcodes) + 1))

# 对于CSR格式, 第i行的数据在 X_data[ (indptr[i]+1) : indptr[i+1] ]
# 列索引在 X_indices[ (indptr[i]+1) : indptr[i+1] ]
# 需要用R的1-based索引

# 只提取匹配细胞的X行
cat(sprintf("  提取 %d 个细胞的数据...\n", length(matched_idx)))

# 分批处理以避免内存问题
n_genes <- h5read(H5AD_FILE, "var/_index")  # gene names
if (length(n_genes) == 0) {
  n_genes <- h5read(H5AD_FILE, "var/gene_symbols")
}
n_genes_total <- length(n_genes)
cat(sprintf("  总基因数: %d\n", n_genes_total))

# 构建稀疏矩阵的行
library(Matrix)
rows_list <- list()
total_nonzero <- 0

for (i in seq_along(matched_idx)) {
  idx <- matched_idx[i]  # 0-based index in h5ad
  start <- X_indptr[idx] + 1  # 转到1-based
  end <- X_indptr[idx + 1]

  if (end >= start) {
    row_data <- X_data[(start:end)]
    row_cols <- X_indices[(start:end)] + 1  # 转到1-based gene index
    rows_list[[i]] <- sparseMatrix(
      i = rep(i, length(row_data)),
      j = row_cols,
      x = row_data,
      dims = c(length(matched_idx), n_genes_total)
    )
    total_nonzero <- total_nonzero + length(row_data)
  }

  if (i %% 500 == 0) {
    cat(sprintf("  已处理: %d/%d cells...\n", i, length(matched_idx)))
  }
}

cat(sprintf("  总非零值: %d\n", total_nonzero))

# 合并所有行
cat("  合并稀疏矩阵...\n")
X_subset <- do.call(rbind, rows_list)
rownames(X_subset) <- all_barcodes[matched_idx]
colnames(X_subset) <- n_genes

cat(sprintf("  子集矩阵: %d genes x %d cells\n", ncol(X_subset), nrow(X_subset)))

# 5. 提取基因信息
cat("\n[4] 提取基因和细胞metadata...\n")

gene_names <- n_genes
# 尝试读gene symbols
tryCatch({
  var <- h5read(H5AD_FILE, "var")
  if ("gene_symbols" %in% names(var)) {
    gene_symbols <- var$gene_symbols
    colnames(X_subset) <- gene_symbols
  }
}, error = function(e) cat(sprintf("  无法读取gene symbols: %s\n", e$message)))

# 6. 保存
cat("\n[5] 保存...\n")
# 保存为R稀疏矩阵
saveRDS(X_subset, OUT_RDS)
cat(sprintf("  表达矩阵: %s\n", OUT_RDS))

# 保存匹配的barcodes
writeLines(all_barcodes[matched_idx],
           "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/matched_barcodes.txt")

cat("\n✅ 提取完成!\n")
cat(sprintf("  细胞数: %d\n", nrow(X_subset)))
cat(sprintf("  基因数: %d\n", ncol(X_subset)))
cat(sprintf("  非零值: %d\n", total_nonzero))
