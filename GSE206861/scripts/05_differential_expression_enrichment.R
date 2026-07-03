###############################################################################
# GSE206861 单细胞分析 — 组间差异分析 + 功能富集
# 比较: Mouse Cystic-fluid (CF) vs Sham
# 重点关注: 小胶质细胞炎症激活、神经元损伤
###############################################################################

library(Seurat)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)
library(here)

setwd(here("GSE206861"))

seu_mouse <- readRDS(here("GSE206861/data/processed/seu_mouse_annotated.rds"))

# ============================================================================
# Part A: 组间差异表达 (CF vs Sham)
# ============================================================================

cat("\n========== Part A: 差异表达分析 (CF vs Sham) ==========\n")

# 设置分组
DefaultAssay(seu_mouse) <- "SCT"
Idents(seu_mouse) <- "orig.ident"

# --- A1. 整体差异基因 ---
cat("\n>>> 整体差异基因 (CF vs Sham)...\n")
de_all <- FindMarkers(
  seu_mouse,
  ident.1 = "Mouse_CF",
  ident.2 = "Mouse_Sham",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
) %>%
  rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC))

cat("  DE genes (all):", nrow(de_all), "\n")
head(de_all, 10)

# --- A2. 每类细胞的差异基因 ---
cat("\n>>> 按细胞类型差异基因...\n")

Idents(seu_mouse) <- seu_mouse$cell_type
cell_types <- unique(seu_mouse$cell_type)

de_by_cell <- list()

for (ct in cell_types) {
  cat("  Cell type:", ct, "\n")

  cells_ct <- WhichCells(seu_mouse, idents = ct)
  cf_cells <- intersect(cells_ct, WhichCells(seu_mouse, expression = orig.ident == "Mouse_CF"))
  sham_cells <- intersect(cells_ct, WhichCells(seu_mouse, expression = orig.ident == "Mouse_Sham"))

  # 两边都需要足够细胞
  if (length(cf_cells) < 10 | length(sham_cells) < 10) {
    cat("    Skipped (not enough cells: CF=", length(cf_cells), "Sham=", length(sham_cells), ")\n")
    next
  }

  de <- FindMarkers(
    seu_mouse,
    ident.1 = cf_cells,
    ident.2 = sham_cells,
    min.pct = 0.1,
    logfc.threshold = 0.25
  ) %>%
    rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC))

  de$cell_type <- ct
  de_by_cell[[ct]] <- de
  cat("    DE genes:", nrow(de), "\n")
}

# 合并所有细胞类型差异结果
de_all_cells <- bind_rows(de_by_cell)
write.csv(de_all_cells, here("GSE206861/results/DE_by_celltype.csv"), row.names = FALSE)

# --- A3. 火山图: 小胶质细胞 ---
if ("Microglia" %in% names(de_by_cell) || any(grepl("icroglia", names(de_by_cell), ignore.case = TRUE))) {
  mg <- names(de_by_cell)[grep("icroglia", names(de_by_cell), ignore.case = TRUE)]
  if (length(mg) > 0) {
    de_mg <- de_by_cell[[mg[1]]]

    # 标注关键基因
    highlight_genes <- c("Cd68", "Cd74", "Il1b", "Tnf", "Apoe", "App", "C1qa",
                         "C1qb", "C1qc", "Trem2", "Tyrobp", "Itgam", "Nlrp3",
                         "Cxcl10", "Ccl2")

    de_mg <- de_mg %>%
      mutate(
        significance = case_when(
          p_val_adj < 0.05 & avg_log2FC > 0.5  ~ "Up (sig)",
          p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down (sig)",
          TRUE ~ "NS"
        ),
        label = ifelse(gene %in% highlight_genes & p_val_adj < 0.05, gene, "")
      )

    p_volcano <- ggplot(de_mg, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                    color = significance, label = label)) +
      geom_point(alpha = 0.6, size = 0.8) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      ggrepel::geom_text_repel(max.overlaps = 30, size = 3, na.rm = TRUE) +
      scale_color_manual(values = c("Up (sig)" = "#E64B35", "Down (sig)" = "#4DBBD5", "NS" = "grey80")) +
      labs(title = "Microglia: CF vs Sham", x = "log2 FC", y = "-log10 adj.P") +
      theme_minimal()

    pdf(here("GSE206861/results/microglia_volcano.pdf"), width = 8, height = 7)
    print(p_volcano)
    dev.off()
  }
}

# --- A4. 热图: 文献关键基因 ---
key_genes <- c("Cd68", "Cd74", "Il1b", "Tnf", "Apoe", "App", "C1qa",
               "Trem2", "Tyrobp", "Itgam", "Nlrp3", "Cxcl10", "Ccl2",
               "Npy", "Fgfr2", "Sst", "Pcsk1n", "Rnpc3", "Agrp")

key_available <- intersect(key_genes, rownames(seu_mouse))

if (length(key_available) > 0) {
  pdf(here("GSE206861/results/key_genes_heatmap.pdf"), width = 14, height = 8)
  DoHeatmap(
    subset(seu_mouse, features = key_available),
    features = key_available,
    group.by = "orig.ident",
    assay = "SCT", slot = "scale.data"
  ) + ggtitle("Key Literature Genes: CF vs Sham")
  dev.off()
}

# ============================================================================
# Part B: 功能富集分析
# ============================================================================

cat("\n========== Part B: 功能富集分析 ==========\n")

# --- B1. 准备基因列表 ---
# 对小胶质细胞差异基因做富集
mg_name <- if (exists("mg") && length(mg) > 0) mg[1] else NULL

if (!is.null(mg_name) && nrow(de_by_cell[[mg_name]]) > 0) {
  de_mg <- de_by_cell[[mg_name]]

  # 上调基因 (CF vs Sham)
  up_genes <- de_mg %>%
    filter(avg_log2FC > 0.5, p_val_adj < 0.05) %>%
    pull(gene)

  # 下调基因
  down_genes <- de_mg %>%
    filter(avg_log2FC < -0.5, p_val_adj < 0.05) %>%
    pull(gene)

  # 背景基因
  bg_genes <- de_mg$gene

  cat("  Microglia UP genes:", length(up_genes), "\n")
  cat("  Microglia DOWN genes:", length(down_genes), "\n")
} else {
  # 降级：用整体差异基因
  cat("  No microglia-specific DE. Using all DE genes.\n")
  up_genes <- de_all %>% filter(avg_log2FC > 0.5, p_val_adj < 0.05) %>% pull(gene)
  down_genes <- de_all %>% filter(avg_log2FC < -0.5, p_val_adj < 0.05) %>% pull(gene)
  bg_genes <- de_all$gene
}

# --- B2. GO 富集 ---
run_go <- function(genes, bg, ont = "BP", title = "") {
  if (length(genes) < 5) {
    cat("  Not enough genes for GO enrichment (", ont, "): ", length(genes), "\n")
    return(NULL)
  }

  ego <- enrichGO(
    gene          = genes,
    universe      = bg,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )

  if (!is.null(ego) && nrow(ego) > 0) {
    cat("  GO", ont, "terms:", nrow(ego), "\n")

    # 绘图
    p_dot <- dotplot(ego, showCategory = 15) + ggtitle(title)
    print(p_dot)
    return(ego)
  } else {
    cat("  No significant GO", ont, "terms.\n")
    return(NULL)
  }
}

cat("\n>>> GO BP - Microglia UP...\n")
ego_up_bp <- run_go(up_genes, bg_genes, "BP", "GO BP: Microglia UP (CF vs Sham)")

cat("\n>>> GO BP - Microglia DOWN...\n")
ego_down_bp <- run_go(down_genes, bg_genes, "BP", "GO BP: Microglia DOWN (CF vs Sham)")

# --- B3. KEGG 富集 ---
run_kegg <- function(genes, bg, title = "") {
  if (length(genes) < 5) return(NULL)

  # 转换为 ENTREZ ID
  gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db, drop = TRUE)
  bg_entrez <- bitr(bg, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db, drop = TRUE)

  if (nrow(gene_entrez) < 5) {
    cat("  Not enough genes after ID mapping.\n")
    return(NULL)
  }

  ekegg <- enrichKEGG(
    gene          = gene_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    organism      = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )

  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    cat("  KEGG pathways:", nrow(ekegg), "\n")
    p_dot <- dotplot(ekegg, showCategory = 15) + ggtitle(title)
    print(p_dot)
    return(ekegg)
  } else {
    cat("  No significant KEGG pathways.\n")
    return(NULL)
  }
}

cat("\n>>> KEGG - Microglia UP...\n")
ekegg_up <- run_kegg(up_genes, bg_genes, "KEGG: Microglia UP (CF vs Sham)")

# --- B4. GSEA — 所有基因排序 ---
cat("\n>>> GSEA...\n")

if (!is.null(mg_name) && nrow(de_by_cell[[mg_name]]) > 0) {
  de_mg <- de_by_cell[[mg_name]]

  # 计算排位指标: -log10(p) * sign(log2FC)
  de_mg <- de_mg %>%
    filter(!is.na(avg_log2FC), !is.na(p_val)) %>%
    mutate(rank_metric = -log10(p_val + 1e-300) * sign(avg_log2FC))

  rank_list <- setNames(de_mg$rank_metric, de_mg$gene)
  rank_list <- sort(rank_list, decreasing = TRUE)

  # 获取 Hallmark 基因集
  hallmark <- msigdbr(species = "Mus musculus", category = "H")
  hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

  gsea_res <- GSEA(
    geneList = rank_list,
    TERM2GENE = hallmark[, c("gs_name", "gene_symbol")],
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500
  )

  if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
    cat("  GSEA significant gene sets:", nrow(gsea_res), "\n")
    print(head(gsea_res[, c("Description", "NES", "p.adjust")], 10))

    pdf(here("GSE206861/results/gsea_hallmark.pdf"), width = 14, height = 10)
    print(dotplot(gsea_res, showCategory = 20) + ggtitle("GSEA: Hallmark (CF vs Sham, Microglia)"))
    dev.off()
  }
}

# --- B5. 保存全部富集结果 ---
dir.create(here("GSE206861/results/enrichment"), recursive = TRUE, showWarnings = FALSE)

if (exists("ego_up_bp") && !is.null(ego_up_bp)) {
  write.csv(as.data.frame(ego_up_bp), here("GSE206861/results/enrichment/GO_BP_up.csv"))
}
if (exists("ego_down_bp") && !is.null(ego_down_bp)) {
  write.csv(as.data.frame(ego_down_bp), here("GSE206861/results/enrichment/GO_BP_down.csv"))
}
if (exists("ekegg_up") && !is.null(ekegg_up)) {
  write.csv(as.data.frame(ekegg_up), here("GSE206861/results/enrichment/KEGG_up.csv"))
}
if (exists("gsea_res") && !is.null(gsea_res)) {
  write.csv(as.data.frame(gsea_res), here("GSE206861/results/enrichment/GSEA_hallmark.csv"))
}

# ============================================================================
# 完成
# ============================================================================
cat("\n✓ 差异分析与富集分析完成！\n")
cat("  结果: results/DE_by_celltype.csv\n")
cat("  结果: results/enrichment/\n")
cat("  接下来运行: 06_advanced_analysis.R\n")
