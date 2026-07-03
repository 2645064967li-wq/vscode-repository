# =============================================================================
# 下丘脑 Microglia 衰老差异分析 (v2 - 使用预提取的小h5ad)
# Jin et al. 2025 Nature
# =============================================================================
set.seed(42)

# ---- 0. 环境 ----
options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_needed <- function(pkgs, type = "cran") {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("安装 %s...\n", pkg))
      if (type == "bioc") {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else if (type == "github") {
        if (!requireNamespace("remotes", quietly = TRUE))
          install.packages("remotes")
        remotes::install_github(pkg, upgrade = "never")
      } else {
        install.packages(pkg)
      }
    }
  }
}

cat("\n============================================\n")
cat(" 下丘脑Microglia 衰老差异分析\n")
cat("============================================\n\n")

# 安装包
install_if_needed(c("Seurat", "dplyr", "ggplot2", "patchwork", "data.table",
                     "Matrix", "cowplot"), "cran")
install_if_needed("mojaveazure/seurat-disk", "github")
install_if_needed("EnhancedVolcano", "bioc")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
})

cat("✓ 包加载完成\n\n")

# ---- 1. 路径 ----
H5AD_SMALL <- "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_microglia_expression.h5ad"
RESULT_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/results"
DATA_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/data"
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 2. 读取子集化后的表达矩阵 ----
cat("[Step 1] 读取下丘脑Microglia表达矩阵...\n")
stopifnot(file.exists(H5AD_SMALL))

# 转换并加载
h5seurat_file <- sub("\\.h5ad$", ".h5seurat", H5AD_SMALL)
if (!file.exists(h5seurat_file)) {
  cat("  转换 h5ad -> h5Seurat...\n")
  Convert(H5AD_SMALL, dest = h5seurat_file, overwrite = TRUE)
}
hypo_micro <- LoadH5Seurat(h5seurat_file)

cat(sprintf("  细胞数: %d\n", ncol(hypo_micro)))
cat(sprintf("  基因数: %d\n", nrow(hypo_micro)))

# ---- 3. 年龄/性别检查 ----
cat("\n[Step 2] 年龄分组...\n")

if ("age_group" %in% colnames(hypo_micro@meta.data)) {
  cat("  年龄分布:\n")
  print(table(hypo_micro$age_group))
  Idents(hypo_micro) <- "age_group"
} else if ("donor_age_category" %in% colnames(hypo_micro@meta.data)) {
  hypo_micro$age_group <- ifelse(hypo_micro$donor_age_category == "aged",
                                  "Aged(18m)", "Adult(2m)")
  Idents(hypo_micro) <- "age_group"
  print(table(hypo_micro$age_group))
} else {
  cat("  ⚠️ 无年龄注释，正在从cell_metadata添加...\n")
  cell_meta <- fread(file.path(DATA_DIR, "single_cell/metadata/cell_metadata.csv"))
  cell_meta <- cell_meta[cell_meta$cell_label %in% colnames(hypo_micro), ]
  hypo_micro$age_group <- ifelse(
    cell_meta$donor_age_category[match(colnames(hypo_micro), cell_meta$cell_label)] == "aged",
    "Aged(18m)", "Adult(2m)")
  Idents(hypo_micro) <- "age_group"
  print(table(hypo_micro$age_group))
}

# ---- 4. 标准化与降维 ----
cat("\n[Step 3] 标准化、高变基因、PCA、UMAP...\n")

# 数据已log2标准化，但仍运行标准流程
hypo_micro <- NormalizeData(hypo_micro, normalization.method = "LogNormalize",
                            scale.factor = 10000)
hypo_micro <- FindVariableFeatures(hypo_micro, selection.method = "vst",
                                    nfeatures = 2000)
cat(sprintf("  高变基因: %d\n", length(VariableFeatures(hypo_micro))))

hypo_micro <- ScaleData(hypo_micro, features = rownames(hypo_micro))
hypo_micro <- RunPCA(hypo_micro, features = VariableFeatures(hypo_micro), npcs = 30)
hypo_micro <- RunUMAP(hypo_micro, dims = 1:15)

# ---- 5. 差异表达: Aged vs Adult ----
cat("\n[Step 4] 差异表达分析...\n")

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

n_sig <- sum(de_results$p_val_adj < 0.05, na.rm = TRUE)
n_up <- sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC > 0, na.rm = TRUE)
n_down <- sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC < 0, na.rm = TRUE)

cat(sprintf("  差异基因 (adj.p<0.05): %d\n", n_sig))
cat(sprintf("    上调(老年): %d\n", n_up))
cat(sprintf("    下调(老年): %d\n", n_down))

# 保存
fwrite(de_results, file.path(RESULT_DIR, "hypo_microglia_aging_DEGs.csv"))
cat(sprintf("  ✓ 已保存: %s\n", file.path(RESULT_DIR, "hypo_microglia_aging_DEGs.csv")))

# Top基因展示
cat("\n  ═══ Top 15 上调 (老年高表达) ═══\n")
top_up <- head(subset(de_results, avg_log2FC > 0), 15)
print(top_up[, c("gene", "avg_log2FC", "p_val_adj")], row.names = FALSE)

cat("\n  ═══ Top 15 下调 (老年低表达) ═══\n")
top_down <- de_results[order(de_results$avg_log2FC), ]
top_down <- head(top_down, 15)
print(top_down[, c("gene", "avg_log2FC", "p_val_adj")], row.names = FALSE)

# ---- 6. AD相关基因 ----
cat("\n[Step 5] AD相关基因检查...\n")

ad_genes <- c(
  "Apoe", "Trem2", "Tyrobp", "Cd33", "Clu", "Bin1", "Picalm",
  "Abca7", "Sorl1", "App", "Bace1", "Psen1", "Psen2", "Adam10",
  "Cd68", "Cd74", "Aif1", "Itgam", "Cx3cr1", "Tmem119", "P2ry12",
  "Cst7", "Lpl", "Spp1", "Clec7a",
  "Il1b", "Il6", "Tnf", "Ccl2", "Ccl3", "C1qa", "C1qb", "C1qc",
  "Mapt", "Snca", "Lrrk2", "Tardbp"
)

ad_in_de <- de_results[de_results$gene %in% ad_genes, ]
if (nrow(ad_in_de) > 0) {
  cat(sprintf("  命中 %d/%d AD相关基因:\n", nrow(ad_in_de), length(ad_genes)))
  print(ad_in_de[order(ad_in_de$avg_log2FC, decreasing = TRUE),
                 c("gene", "avg_log2FC", "p_val_adj")], row.names = FALSE)
} else {
  cat("  ⚠️ 基因命名格式可能不同，请检查基因符号格式\n")
  cat("  示例基因名 (前10个):\n")
  print(head(rownames(hypo_micro), 10))
}

# ---- 7. GO富集 (简化版, 使用top DEGs) ----
cat("\n[Step 6] GO富集分析...\n")

if (n_sig >= 10) {
  sig_genes <- rownames(subset(de_results, p_val_adj < 0.05))
  cat(sprintf("  显著性基因: %d (用于富集)\n", length(sig_genes)))

  # 如果用clusterProfiler
  install_if_needed(c("clusterProfiler", "org.Mm.eg.db", "enrichplot"), "bioc")

  if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    suppressPackageStartupMessages({
      library(clusterProfiler)
      library(org.Mm.eg.db)
      library(enrichplot)
    })

    # 转换gene symbol -> ENTREZ ID
    sig_genes_clean <- gsub("-", ".", sig_genes)  # Seurat会转换_为-

    # 尝试GO BP富集 (可能需要symbol转换)
    tryCatch({
      go_results <- enrichGO(
        gene = sig_genes[1:min(500, length(sig_genes))],
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )

      if (!is.null(go_results) && nrow(go_results) > 0) {
        go_simple <- go_results@result[, c("Description", "p.adjust", "Count")]
        go_simple <- go_simple[order(go_simple$p.adjust), ]
        cat(sprintf("  富集GO terms: %d\n", nrow(go_simple)))
        cat("\n  Top 20 GO terms:\n")
        print(head(go_simple, 20), row.names = FALSE)

        fwrite(go_simple, file.path(RESULT_DIR, "hypo_micro_GO_enrichment.csv"))
      } else {
        cat("  GO富集无显著结果\n")
      }
    }, error = function(e) {
      cat(sprintf("  GO富集失败: %s\n", e$message))
    })
  }
}

# ---- 8. 可视化 ----
cat("\n[Step 7] 可视化...\n")

# UMAP - 年龄
p1 <- DimPlot(hypo_micro, group.by = "age_group",
              cols = c("Adult(2m)" = "#4DBBD5", "Aged(18m)" = "#E64B35"),
              pt.size = 1.5) +
  ggtitle("Hypothalamic Microglia: Age Groups") +
  theme_minimal(base_size = 12)

# UMAP - cluster
if ("cluster_alias" %in% colnames(hypo_micro@meta.data)) {
  hypo_micro$cluster_alias <- as.factor(hypo_micro$cluster_alias)
  p2 <- DimPlot(hypo_micro, group.by = "cluster_alias", pt.size = 1.5) +
    ggtitle("Hypothalamic Microglia: Subclusters") +
    theme_minimal(base_size = 12)
}

# 火山图
de_plot <- de_results[!is.na(de_results$p_val_adj), ]
de_plot$signif <- "Not Sig"
de_plot$signif[de_plot$p_val_adj < 0.05 & de_plot$avg_log2FC > 0.5] <- "Up in Aged"
de_plot$signif[de_plot$p_val_adj < 0.05 & de_plot$avg_log2FC < -0.5] <- "Down in Aged"

# 标记显著AD基因
de_plot$label <- ""
ad_high <- de_plot$gene %in% ad_genes & de_plot$p_val_adj < 0.05
de_plot$label[ad_high] <- de_plot$gene[ad_high]

p3 <- ggplot(de_plot, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300),
                           color = signif)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = c("Up in Aged" = "#E64B35",
                                 "Down in Aged" = "#4DBBD5",
                                 "Not Sig" = "grey85")) +
  ggrepel::geom_text_repel(aes(label = label), size = 3.5,
                            max.overlaps = 30, color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
  ggtitle("Hypothalamic Microglia: Aged vs Adult") +
  xlab("log2 Fold Change (Aged/Adult)") + ylab("-log10 adjusted p-value") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

# Top DEGs 热图
top_degs <- unique(c(
  head(subset(de_results, avg_log2FC > 0 & p_val_adj < 0.05)$gene, 15),
  head(subset(de_results, avg_log2FC < 0 & p_val_adj < 0.05)$gene, 15)
))

if (length(top_degs) >= 4) {
  p4 <- DoHeatmap(hypo_micro, features = top_degs, group.by = "age_group",
                   size = 3.5, disp.min = -2, disp.max = 2) +
    ggtitle("Top DEGs: Hypothalamic Microglia Aging") +
    theme(axis.text.y = element_text(size = 7))
} else {
  p4 <- NULL
}

# 保存图
ggsave(file.path(RESULT_DIR, "01_umap_age.png"), p1, width = 7, height = 6, dpi = 150)
if (exists("p2")) ggsave(file.path(RESULT_DIR, "02_umap_cluster.png"), p2, width = 8, height = 6, dpi = 150)
ggsave(file.path(RESULT_DIR, "03_volcano.png"), p3, width = 10, height = 8, dpi = 150)
if (!is.null(p4)) ggsave(file.path(RESULT_DIR, "04_heatmap.png"), p4, width = 12, height = 8, dpi = 150)
cat("  ✓ 图表已保存\n")

# ---- 9. 保存对象 ----
cat("\n[Step 8] 保存...\n")
saveRDS(hypo_micro, file.path(DATA_DIR, "hypothalamus/hypo_microglia_seurat.rds"))
cat("  ✓ Seurat RDS已保存\n")

# ---- 10. 总结 ----
cat("\n============================================\n")
cat(" ✅ 分析完成!\n")
cat("============================================\n")
cat(sprintf(" 细胞数: %d\n", ncol(hypo_micro)))
cat(sprintf(" DEGs (adj.p<0.05): %d (↑%d ↓%d)\n", n_sig, n_up, n_down))
if (exists("go_results") && !is.null(go_results)) {
  cat(sprintf(" GO terms: %d\n", nrow(go_results)))
}
cat(sprintf("\n 结果: %s/\n", RESULT_DIR))
cat("  - hypo_microglia_aging_DEGs.csv\n")
cat("  - hypo_micro_GO_enrichment.csv\n")
cat("  - hypo_microglia_seurat.rds\n")
