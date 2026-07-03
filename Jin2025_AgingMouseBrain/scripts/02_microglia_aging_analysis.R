# =============================================================================
# 下丘脑 Microglia 衰老差异分析
# Jin et al. 2025 Nature - Brain-wide ageing signatures in mice
# =============================================================================
# 分析流程:
#   1. 安装/加载依赖包
#   2. 读取表达矩阵 (h5ad)
#   3. 筛选下丘脑Microglia
#   4. 成年 vs 老年 差异表达分析
#   5. AD相关基因探索
# =============================================================================

# =============================================================================
# 0. 环境配置 ----
# =============================================================================
set.seed(42)

# 镜像设置
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

# 包安装函数
install_if_needed <- function(pkgs, type = "cran") {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("安装 %s (%s)...\n", pkg, type))
      if (type == "cran") {
        install.packages(pkg, quiet = FALSE)
      } else if (type == "bioc") {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager", quiet = TRUE)
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else if (type == "github") {
        if (!requireNamespace("remotes", quietly = TRUE))
          install.packages("remotes", quiet = TRUE)
        remotes::install_github(pkg, upgrade = "never")
      }
    }
  }
}

cat("\n============================================\n")
cat(" 下丘脑 Microglia 衰老差异分析\n")
cat(" Jin et al. 2025 Nature\n")
cat("============================================\n\n")

# Step 0a: CRAN基础包
cran_pkgs <- c("Seurat", "Matrix", "dplyr", "ggplot2", "patchwork",
               "cowplot", "remotes", "BiocManager", "data.table")
install_if_needed(cran_pkgs, "cran")

# Step 0b: SeuratDisk (GitHub only, CRAN version removed)
install_if_needed("mojaveazure/seurat-disk", "github")

# Step 0c: EnhancedVolcano (Bioconductor)
install_if_needed("EnhancedVolcano", "bioc")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

cat("✓ 包加载完成\n\n")

# =============================================================================
# 1. 路径设置 ----
# =============================================================================
H5AD_FILE <- "d:/decrepitude mouse hypothamulas/Zeng-Aging-Mouse-10Xv3-log2.h5ad"
BARCODE_FILE <- "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_microglia_barcodes.txt"
RESULT_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/results"
DATA_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/data"

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("表达矩阵: %s\n", H5AD_FILE))
cat(sprintf("Microglia barcodes: %s\n", BARCODE_FILE))

# 检查文件
stopifnot(file.exists(H5AD_FILE))
stopifnot(file.exists(BARCODE_FILE))

# 读入下丘脑microglia barcodes
hypo_micro_barcodes <- readLines(BARCODE_FILE)
cat(sprintf("下丘脑Microglia barcodes: %d\n", length(hypo_micro_barcodes)))

# =============================================================================
# 2. 读取 h5ad 表达矩阵 ----
# =============================================================================
cat("\n[Step 1] 读取 h5ad 文件...\n")
cat("  (13GB文件，加载需要5-15分钟，请耐心等待)\n")

# 方法A: 直接读取h5ad (SeuratDisk需要)
# 先转换为h5Seurat格式
h5seurat_file <- file.path(DATA_DIR, "single_cell/expression/Zeng-Aging-Mouse-10Xv3.h5seurat")

if (!file.exists(h5seurat_file)) {
  cat("  转换 h5ad -> h5Seurat (首次需要)...\n")
  Convert(H5AD_FILE, dest = h5seurat_file, overwrite = FALSE)
  cat("  ✓ 转换完成\n")
}

cat("  加载 h5Seurat...\n")
seurat_obj <- LoadH5Seurat(h5seurat_file)

cat(sprintf("  总细胞: %s\n", format(ncol(seurat_obj), big.mark = ",")))
cat(sprintf("  总基因: %s\n", format(nrow(seurat_obj), big.mark = ",")))

# =============================================================================
# 3. 筛选下丘脑Microglia ----
# =============================================================================
cat("\n[Step 2] 筛选下丘脑Microglia...\n")

# 匹配barcodes (注意: h5ad中的barcode格式可能与txt中的略有不同)
available_barcodes <- colnames(seurat_obj)
matched_barcodes <- intersect(hypo_micro_barcodes, available_barcodes)

cat(sprintf("  匹配到: %d / %d cells\n",
            length(matched_barcodes), length(hypo_micro_barcodes)))

if (length(matched_barcodes) == 0) {
  cat("\n  ⚠️ barcode格式不匹配! 尝试模糊匹配...\n")
  cat("  前5个available barcodes:\n")
  print(head(available_barcodes, 5))
  cat("  前5个target barcodes:\n")
  print(head(hypo_micro_barcodes, 5))

  # 可能前缀/后缀不同，尝试匹配第一个
  barcode_starts <- gsub("-.*$", "", hypo_micro_barcodes[1:3])
  cat(sprintf("  搜索: %s\n", barcode_starts[1]))
  matches <- grep(barcode_starts[1], available_barcodes, value = TRUE)
  cat(sprintf("  找到: %s\n", paste(head(matches, 3), collapse = ", ")))
}

# 子集化
if (length(matched_barcodes) > 0) {
  hypo_micro <- subset(seurat_obj, cells = matched_barcodes)
  cat(sprintf("  ✓ 下丘脑Microglia子集: %d cells\n", ncol(hypo_micro)))
}

# 如果直接匹配失败，使用cluster信息筛选
if (length(matched_barcodes) == 0) {
  cat("\n  使用cluster信息筛选...\n")
  # 检查是否有cluster信息在metadata中
  if ("cluster_alias" %in% colnames(seurat_obj@meta.data)) {
    # Microglia clusters: 840-844, 下丘脑用解剖分区
    if ("anatomical_division_label" %in% colnames(seurat_obj@meta.data)) {
      hypo_micro <- subset(seurat_obj,
        cluster_alias >= 840 & cluster_alias <= 844 &
        anatomical_division_label == "HY - HY")
    } else {
      hypo_micro <- subset(seurat_obj,
        cluster_alias >= 840 & cluster_alias <= 844)
    }
    cat(sprintf("  ✓ 通过cluster筛选: %d cells\n", ncol(hypo_micro)))
  }
}

# 释放全脑数据，节省内存
rm(seurat_obj)
gc()

# =============================================================================
# 4. 添加年龄和性别注释 ----
# =============================================================================
cat("\n[Step 3] 整理年龄/性别注释...\n")

# 从cell_metadata中添加年龄信息
cell_meta_file <- file.path(DATA_DIR, "single_cell/metadata/cell_metadata.csv")

if (file.exists(cell_meta_file)) {
  cat("  加载cell_metadata...\n")
  cell_meta <- data.table::fread(cell_meta_file, showProgress = FALSE)

  # 筛选对应的cells
  target_barcodes <- colnames(hypo_micro)
  cell_meta_sub <- cell_meta[cell_meta$cell_label %in% target_barcodes, ]
  cat(sprintf("  匹配metadata: %d / %d\n", nrow(cell_meta_sub), ncol(hypo_micro)))

  # 添加到Seurat对象
  # 用donor信息获取年龄
  donor_file <- file.path(DATA_DIR, "single_cell/metadata/donor.csv")
  donor_info <- data.table::fread(donor_file)

  # 简化annotations
  cell_meta_sub$age_group <- ifelse(
    cell_meta_sub$donor_age_category == "adult", "Adult(2m)", "Aged(18m)"
  )

  # 添加到meta.data
  rownames(cell_meta_sub) <- cell_meta_sub$cell_label
  common_cells <- intersect(colnames(hypo_micro), cell_meta_sub$cell_label)

  if (length(common_cells) > 0) {
    hypo_micro$age_group <- cell_meta_sub[common_cells, ]$age_group
    hypo_micro$donor_sex <- cell_meta_sub[common_cells, ]$donor_sex
    hypo_micro$cluster_alias <- cell_meta_sub[common_cells, ]$cluster_alias
    cat(sprintf("  ✓ 添加了年龄/性别注释: %d cells\n", length(common_cells)))

    cat("\n  年龄分布:\n")
    print(table(hypo_micro$age_group))
    cat("\n  性别x年龄:\n")
    print(table(hypo_micro$donor_sex, hypo_micro$age_group))
  }
}

# =============================================================================
# 5. 数据标准化和降维 ----
# =============================================================================
cat("\n[Step 4] 数据标准化和降维...\n")

# 即使h5ad中已有log2值，也运行Seurat的标准流程
hypo_micro <- NormalizeData(hypo_micro, normalization.method = "LogNormalize",
                            scale.factor = 10000)

# 找高变基因
hypo_micro <- FindVariableFeatures(hypo_micro, selection.method = "vst",
                                    nfeatures = 2000)
cat(sprintf("  高变基因: %d\n", length(VariableFeatures(hypo_micro))))

# Scale
hypo_micro <- ScaleData(hypo_micro, features = rownames(hypo_micro))

# PCA
hypo_micro <- RunPCA(hypo_micro, features = VariableFeatures(hypo_micro), npcs = 30)

# UMAP
hypo_micro <- RunUMAP(hypo_micro, dims = 1:15)

# =============================================================================
# 6. 成年 vs 老年 差异表达分析 ----
# =============================================================================
cat("\n[Step 5] 差异表达分析: Adult(2m) vs Aged(18m)...\n")

# 设置细胞分组
Idents(hypo_micro) <- "age_group"

# 找差异基因 (老年 vs 成年)
de_results <- FindMarkers(
  hypo_micro,
  ident.1 = "Aged(18m)",
  ident.2 = "Adult(2m)",
  min.pct = 0.1,
  logfc.threshold = 0.1,
  test.use = "wilcox"
)

de_results$gene <- rownames(de_results)
de_results <- de_results[order(de_results$avg_log2FC, decreasing = TRUE), ]

cat(sprintf("  差异基因总数: %d (p_val_adj < 0.05)\n",
            sum(de_results$p_val_adj < 0.05, na.rm = TRUE)))
cat(sprintf("  上调(老年): %d\n",
            sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC > 0, na.rm = TRUE)))
cat(sprintf("  下调(老年): %d\n",
            sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC < 0, na.rm = TRUE)))

# 保存结果
de_file <- file.path(RESULT_DIR, "hypo_microglia_aging_DEGs.csv")
write.csv(de_results, de_file, row.names = FALSE)
cat(sprintf("  ✓ 保存: %s\n", de_file))

# Top DEGs
cat("\n  === Top 20 上调基因 (老年) ===\n")
top_up <- head(subset(de_results, avg_log2FC > 0), 20)
print(top_up[, c("gene", "avg_log2FC", "p_val_adj")])

cat("\n  === Top 20 下调基因 (老年) ===\n")
top_down <- head(subset(de_results, avg_log2FC < 0), 20)
top_down <- top_down[order(top_down$avg_log2FC), ]
print(top_down[, c("gene", "avg_log2FC", "p_val_adj")])

# =============================================================================
# 7. AD相关基因分析 ----
# =============================================================================
cat("\n[Step 6] AD相关基因分析...\n")

# 已知AD风险基因
ad_genes <- c(
  # GWAS风险基因
  "Apoe", "Trem2", "Tyrobp", "Cd33", "Clu", "Bin1", "Picalm",
  "Abca7", "Sorl1", "Bace1", "Psen1", "Psen2", "App", "Adam10",
  # 小胶质细胞激活
  "Cd68", "Cd74", "Aif1", "Itgam", "Cx3cr1", "Tmem119", "P2ry12",
  "Cst7", "Lpl", "Spp1", "Clec7a",
  # 炎症因子
  "Il1b", "Il6", "Tnf", "Ccl2", "Ccl3", "C1qa", "C1qb", "C1qc",
  # 衰老/神经退行
  "Mapt", "Snca", "Lrrk2", "Tardbp"
)

# 检查这些基因在差异表达中的情况
ad_in_de <- de_results[de_results$gene %in% ad_genes, ]
if (nrow(ad_in_de) > 0) {
  cat(sprintf("  AD基因在DEGs中: %d / %d\n", nrow(ad_in_de), length(ad_genes)))
  cat("\n  AD基因差异表达:\n")
  print(ad_in_de[, c("gene", "avg_log2FC", "p_val_adj")])
} else {
  cat("  ⚠️ 基因名可能使用不同命名格式\n")
}

# =============================================================================
# 8. 可视化 ----
# =============================================================================
cat("\n[Step 7] 生成可视化图表...\n")

# 8.1 UMAP按年龄分组
p1 <- DimPlot(hypo_micro, group.by = "age_group",
              cols = c("Adult(2m)" = "#4DBBD5", "Aged(18m)" = "#E64B35"),
              pt.size = 1) +
  ggtitle("下丘脑Microglia: 年龄分组") +
  theme_minimal()

# 8.2 UMAP按microglia亚群
if ("cluster_alias" %in% colnames(hypo_micro@meta.data)) {
  p2 <- DimPlot(hypo_micro, group.by = "cluster_alias", pt.size = 1) +
    ggtitle("下丘脑Microglia: 亚群") +
    theme_minimal()
}

# 8.3 Top DEGs热图
top_degs <- unique(c(
  head(subset(de_results, avg_log2FC > 0 & p_val_adj < 0.05)$gene, 15),
  head(subset(de_results, avg_log2FC < 0 & p_val_adj < 0.05)$gene, 15)
))

if (length(top_degs) > 2) {
  p3 <- DoHeatmap(hypo_micro, features = top_degs, group.by = "age_group",
                   size = 3) +
    ggtitle("Top 差异基因 (Aged vs Adult)") +
    theme(axis.text.y = element_text(size = 6))
} else {
  p3 <- NULL
  cat("  ⚠️ 差异基因不足，跳过热图\n")
}

# 8.4 火山图 (基础版)
de_plot <- de_results[!is.na(de_results$p_val_adj), ]
de_plot$significance <- "Not Significant"
de_plot$significance[de_plot$p_val_adj < 0.05 & de_plot$avg_log2FC > 0.5] <- "Up in Aged"
de_plot$significance[de_plot$p_val_adj < 0.05 & de_plot$avg_log2FC < -0.5] <- "Down in Aged"

# 标记AD基因
de_plot$label <- ifelse(de_plot$gene %in% ad_genes, de_plot$gene, "")

p4 <- ggplot(de_plot, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300),
                           color = significance)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("Up in Aged" = "#E64B35",
                                 "Down in Aged" = "#4DBBD5",
                                 "Not Significant" = "grey80")) +
  geom_text(aes(label = label), size = 3, vjust = -1, color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
  ggtitle("下丘脑Microglia: Aged vs Adult 火山图") +
  xlab("log2 Fold Change (Aged/Adult)") +
  ylab("-log10 adjusted p-value") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 保存图片
ggsave(file.path(RESULT_DIR, "01_umap_age.png"), p1, width = 8, height = 6, dpi = 150)
if (exists("p2")) ggsave(file.path(RESULT_DIR, "02_umap_cluster.png"), p2, width = 8, height = 6, dpi = 150)
if (!is.null(p3)) ggsave(file.path(RESULT_DIR, "03_heatmap.png"), p3, width = 12, height = 8, dpi = 150)
ggsave(file.path(RESULT_DIR, "04_volcano.png"), p4, width = 10, height = 8, dpi = 150)

cat("  ✓ 图表已保存到 results/\n")

# =============================================================================
# 9. 保存R对象 ----
# =============================================================================
cat("\n[Step 8] 保存分析对象...\n")

rds_file <- file.path(DATA_DIR, "hypothalamus/hypo_microglia_seurat.rds")
saveRDS(hypo_micro, rds_file)
cat(sprintf("  ✓ Seurat对象: %s\n", rds_file))

# =============================================================================
# 10. 报告总结 ----
# =============================================================================
cat("\n============================================\n")
cat(" 分析完成! 结果摘要\n")
cat("============================================\n")
cat(sprintf(" 下丘脑Microglia: %d cells\n", ncol(hypo_micro)))
cat(sprintf(" 差异基因(adj.p<0.05): %d\n", sum(de_results$p_val_adj < 0.05, na.rm = TRUE)))
cat(sprintf(" 上调(Aged): %d | 下调(Aged): %d\n",
    sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC > 0, na.rm = TRUE),
    sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC < 0, na.rm = TRUE)))
cat(sprintf(" AD风险基因差异表达: %d genes\n", nrow(ad_in_de)))

cat("\n输出文件:\n")
cat(sprintf("  差异基因表: %s\n", de_file))
cat(sprintf("  Seurat对象: %s\n", rds_file))
cat(sprintf("  图表: %s/\n", RESULT_DIR))
cat("    - 01_umap_age.png\n")
cat("    - 02_umap_cluster.png\n")
cat("    - 03_heatmap.png\n")
cat("    - 04_volcano.png\n")

cat("\n下一步建议:\n")
cat("  1. 查看差异基因列表,关注AD相关基因\n")
cat("  2. GO/KEGG富集分析\n")
cat("  3. 比较不同microglia亚群的差异\n")
cat("  4. 与MERFISH空间数据整合\n")
