###############################################################################
# GSE206861 单细胞分析 — 进阶: 拟时序 + 细胞通讯
# 拟时序: Monocle3 推断小胶质细胞激活轨迹
# 细胞通讯: CellChat 分析 CD74-APP 等配体-受体互作
###############################################################################

library(Seurat)
library(tidyverse)
library(patchwork)
library(here)

setwd(here("GSE206861"))

seu_mouse <- readRDS(here("GSE206861/data/processed/seu_mouse_annotated.rds"))

# ============================================================================
# Part A: 拟时序分析 (Monocle3) — 小胶质细胞激活轨迹
# ============================================================================

cat("\n========== Part A: 拟时序分析 (Monocle3) ==========\n")

library(monocle3)
library(SeuratWrappers)

# --- A1. 提取小胶质细胞 ---
mg_idx <- grep("icroglia", unique(seu_mouse$cell_type), ignore.case = TRUE, value = TRUE)

if (length(mg_idx) > 0) {
  cat(">>> 提取小胶质细胞...\n")
  seu_mg <- subset(seu_mouse, cell_type %in% mg_idx)

  # 重新降维 (用更少 PC 捕获连续轨迹)
  seu_mg <- NormalizeData(seu_mg, assay = "RNA") %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(npcs = 30)

  cat("  Microglia cells:", ncol(seu_mg), "\n")

  # --- A2. 转换为 Monocle3 CellDataSet ---
  cds <- as.cell_data_set(seu_mg, assay = "RNA")

  # 转换变量名
  colData(cds)$sample <- seu_mg$orig.ident
  colData(cds)$cluster <- seu_mg$seurat_clusters

  # --- A3. 拟时序: 预处理 + 降维 + 排序 ---
  cat(">>> Monocle3 preprocessing...\n")
  cds <- preprocess_cds(cds, num_dim = 20, method = "PCA")

  cat(">>> UMAP (Monocle3)...\n")
  cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP")

  cat(">>> 学习轨迹图...\n")
  cds <- cluster_cells(cds, resolution = 1e-3)
  cds <- learn_graph(cds)

  cat(">>> 排序伪时间...\n")
  # 选择根节点: 用稳态小胶质细胞 marker 表达最高的 cluster
  # Tmem119, P2ry12 高表达 = 稳态
  cds <- order_cells(cds)

  # --- A4. 可视化 ---
  pdf(here("GSE206861/results/microglia_trajectory.pdf"), width = 18, height = 6)

  p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                   label_branch_points = FALSE, label_roots = FALSE,
                   label_leaves = FALSE, graph_label_size = 2) +
    ggtitle("Microglia: Pseudotime Trajectory")

  p2 <- plot_cells(cds, color_cells_by = "sample", label_cell_groups = FALSE) +
    ggtitle("Microglia: Sample (CF=Mouse_CF, Sham=Mouse_Sham)")

  p3 <- plot_cells(cds, color_cells_by = "cluster", label_cell_groups = TRUE,
                   group_label_size = 3) +
    ggtitle("Microglia: Clusters")

  print(p1 + p2 + p3 + plot_layout(ncol = 3))
  dev.off()

  # --- A5. 沿伪时间的基因表达 ---
  # 炎症 marker
  inflam_genes <- intersect(c("Cd68", "Cd74", "Il1b", "Tnf", "Apoe", "Trem2", "C1qa"),
                            rownames(cds))

  p_time <- NULL
  if (length(inflam_genes) > 0) {
    pdf(here("GSE206861/results/microglia_pseudotime_genes.pdf"), width = 16, height = 10)
    for (g in inflam_genes) {
      p <- plot_cells(cds, genes = g, label_cell_groups = FALSE,
                      show_trajectory_graph = FALSE) + ggtitle(g)
      p_time <- if (is.null(p_time)) p else p_time + p
    }
    print(p_time + plot_layout(ncol = 3))
    dev.off()
  }

  # --- A6. 保存 ---
  saveRDS(cds, here("GSE206861/data/processed/microglia_cds.rds"))
  cat("  Monocle3 轨迹分析完成。\n")

} else {
  cat("  !! 未找到小胶质细胞注释，跳过拟时序分析。\n")
}

# ============================================================================
# Part B: 细胞通讯分析 (CellChat)
# ============================================================================

cat("\n========== Part B: 细胞通讯分析 (CellChat) ==========\n")

library(CellChat)

# --- B1. 提取小鼠数据（RNA assay）---
DefaultAssay(seu_mouse) <- "RNA"

# 需要 counts 数据和 meta 数据
cellchat_input <- GetAssayData(seu_mouse, layer = "data")  # normalized data

# --- B2. 分别创建 CF 和 Sham 的 CellChat 对象 ---
meta <- seu_mouse@meta.data

create_cellchat <- function(seu_obj, sample_name, species = "mouse") {
  cells <- colnames(seu_obj)[seu_obj$orig.ident == sample_name]

  if (length(cells) < 100) {
    cat("  !! Not enough cells for", sample_name, ":", length(cells), "\n")
    return(NULL)
  }

  data_mat <- GetAssayData(seu_obj, layer = "data")[, cells, drop = FALSE]
  meta_sub <- seu_obj@meta.data[cells, ]

  cellchat <- createCellChat(object = data_mat, meta = meta_sub, group.by = "cell_type")

  # 设置配体-受体数据库
  if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
  } else {
    CellChatDB <- CellChatDB.human
  }

  cellchat@DB <- CellChatDB

  # 预处理
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # 通讯概率计算
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)
}

cat("\n>>> CellChat: Cystic-fluid group...\n")
cc_cf <- create_cellchat(seu_mouse, "Mouse_CF", "mouse")

cat("\n>>> CellChat: Sham group...\n")
cc_sham <- create_cellchat(seu_mouse, "Mouse_Sham", "mouse")

# --- B3. 比较分析 ---
if (!is.null(cc_cf) && !is.null(cc_sham)) {
  cat("\n>>> 合并 CellChat 对象...\n")
  cc_list <- list(CF = cc_cf, Sham = cc_sham)
  cc_merged <- mergeCellChat(cc_list, add.names = names(cc_list))

  # --- B4. 可视化 ---
  pdf(here("GSE206861/results/cellchat_interaction_numbers.pdf"), width = 12, height = 6)
  compareInteractions(cc_merged, show.legend = FALSE, group = c(1, 2))
  dev.off()

  pdf(here("GSE206861/results/cellchat_circle.pdf"), width = 16, height = 8)
  par(mfrow = c(1, 2))
  netVisual_circle(cc_cf@net$weight, vertex.weight = as.numeric(table(cc_cf@idents)),
                   title.name = "CF: Interaction strength")
  netVisual_circle(cc_sham@net$weight, vertex.weight = as.numeric(table(cc_sham@idents)),
                   title.name = "Sham: Interaction strength")
  dev.off()

  # --- B5. 热图: 差异信号通路 ---
  pdf(here("GSE206861/results/cellchat_pathway_heatmap.pdf"), width = 12, height = 10)
  rankNet(cc_merged, mode = "comparison", measure = "weight", stacked = FALSE,
          color.use = c("#E64B35", "#4DBBD5"))
  dev.off()

  # --- B6. 关注 CD74-APP 通路 ---
  cat("\n>>> 检测 CD74-APP 互作...\n")

  # 检查 CD74-APP 是否在配体列表中
  search_pair <- function(cc, source = "CD74", target = "APP") {
    lr_pair <- paste(source, target, sep = "_")
    all_interactions <- cc@LR$LRsig$interaction_name
    grep(lr_pair, all_interactions, value = TRUE, ignore.case = TRUE)
  }

  cat("  CF group CD74-APP interactions:\n")
  if (!is.null(cc_cf)) {
    cd74_cf <- search_pair(cc_cf)
    if (length(cd74_cf) > 0) print(cd74_cf) else cat("  (none found with exact name)\n")
  }

  cat("  Sham group CD74-APP interactions:\n")
  if (!is.null(cc_sham)) {
    cd74_sham <- search_pair(cc_sham)
    if (length(cd74_sham) > 0) print(cd74_sham) else cat("  (none found with exact name)\n")
  }

  # --- B7. 特定通路可视化 ---
  # MIF-CD74, APP-CD74 等通路
  key_pathways <- intersect(c("MIF", "APP", "CD74", "CD44", "CXCL", "CCL", "IL1", "TNF"),
                            names(cc_cf@netP$pathways))

  if (length(key_pathways) > 0) {
    pdf(here("GSE206861/results/cellchat_key_pathways.pdf"), width = 14, height = 10)
    for (pw in key_pathways) {
      netVisual_aggregate(cc_cf, signaling = pw, layout = "circle",
                         vertex.label.cex = 0.7)
      title(main = paste("CF:", pw))
    }
    dev.off()
  }

  # --- B8. 保存 ---
  saveRDS(cc_merged, here("GSE206861/data/processed/cellchat_merged.rds"))
}

# ============================================================================
# 完成
# ============================================================================
cat("\n✓ 进阶分析完成！\n")
cat("  保存: data/processed/microglia_cds.rds (Monocle3)\n")
cat("  保存: data/processed/cellchat_merged.rds (CellChat)\n")
cat("  接下来运行: 07_master_run.R 或查看 results/ 中的图\n")
