###############################################################################
# GSE206861 单细胞分析 — 细胞注释
# 自动注释 (SingleR) + 手动 Marker 验证
# Neuron, Microglia, Astrocyte, Oligodendrocyte, Endothelial, Pericyte, etc.
###############################################################################

library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(celldex)
library(here)

setwd(here("GSE206861"))
set.seed(42)

# 读取整合后数据
seu_mouse <- readRDS(here("GSE206861/data/processed/seu_mouse_integrated.rds"))
seu_acp   <- readRDS(here("GSE206861/data/processed/seu_acp_integrated.rds"))

# ============================================================================
# Part A: 小鼠样本注释
# ============================================================================

cat("\n========== 小鼠细胞注释 ==========\n")

# --- A1. 载入参考数据集 ---
cat("\n>>> 加载参考数据...\n")
ref_mouse <- MouseRNAseqData()  # ImmGen 参考，主要免疫细胞 + 一些基质细胞

# 如果需要更全面的脑组织参考，可以用 Allen Brain Atlas
# ref_brain <- celldex::BlueprintEncodeData()  # 包含更多基质细胞

# --- A2. SingleR 自动注释 ---
cat("\n>>> SingleR 注释...\n")
# 从 SCT 中提取 counts 用于 SingleR
DefaultAssay(seu_mouse) <- "SCT"
mouse_counts <- GetAssayData(seu_mouse, layer = "counts")

pred_mouse <- SingleR(
  test   = mouse_counts,
  ref    = ref_mouse,
  labels = ref_mouse$label.main,    # 主标签
  de.method = "classic"
)

seu_mouse$singleR_main <- pred_mouse$labels[match(colnames(seu_mouse), rownames(pred_mouse))]
seu_mouse$singleR_score <- pred_mouse$tuning.scores$first[match(colnames(seu_mouse), rownames(pred_mouse))]

cat("  SingleR 注释结果:\n")
print(table(seu_mouse$singleR_main))

# --- A3. Manual marker 验证 ---
cat("\n>>> Marker 基因表达验证...\n")

# 下丘脑细胞 marker 列表
markers <- list(
  "Neuron"        = c("Rbfox3", "Syt1", "Tubb3", "Snap25", "Syp"),
  "Microglia"     = c("Cx3cr1", "Tmem119", "P2ry12", "Csf1r", "Aif1"),
  "Astrocyte"     = c("Gfap", "Aqp4", "Aldh1l1", "Slc1a3", "Glast"),
  "Oligodendrocyte" = c("Mbp", "Mog", "Olig1", "Olig2", "Plp1"),
  "OPC"           = c("Pdgfra", "Cspg4", "Sox10"),
  "Endothelial"   = c("Cldn5", "Pecam1", "Flt1", "Cdh5"),
  "Pericyte"      = c("Pdgfrb", "Rgs5", "Cspg4", "Anpep"),
  "Ependymal"     = c("Foxj1", "Cfap54", "Dynlrb2"),
  "Tanycyte"      = c("Rax", "Dio2", "Vim"),  # 下丘脑特殊细胞
  "T_cell"        = c("Cd3d", "Cd3e", "Cd3g", "Trbc2"),
  "B_cell"        = c("Cd79a", "Cd79b", "Ms4a1", "Cd19"),
  "NK_cell"       = c("Nkg7", "Klrb1c", "Klrk1")
)

# 检查哪些 marker 存在于数据中
all_markers <- unique(unlist(markers))
markers_available <- all_markers[all_markers %in% rownames(seu_mouse)]
markers_missing <- setdiff(all_markers, markers_available)
cat("  Available markers:", length(markers_available), "/", length(all_markers), "\n")
if (length(markers_missing) > 0) {
  cat("  Missing markers:", paste(markers_missing, collapse = ", "), "\n")
}

# --- A4. DotPlot ---
pdf(here("GSE206861/results/mouse_markers_dotplot.pdf"), width = 20, height = 8)
DotPlot(seu_mouse, features = unique(unlist(markers)), assay = "SCT",
        group.by = "seurat_clusters") +
  RotatedAxis() +
  ggtitle("Mouse: Marker Gene Expression by Cluster")
dev.off()

# --- A5. Feature Plot: 核心细胞类型 ---
core_features <- c("Rbfox3", "Cx3cr1", "Tmem119", "Gfap", "Aqp4", "Mbp",
                   "Cldn5", "Pdgfra", "Pdgfrb", "Cd3d", "Aif1")
core_available <- intersect(core_features, rownames(seu_mouse))

pdf(here("GSE206861/results/mouse_feature_core.pdf"), width = 16, height = 14)
FeaturePlot(seu_mouse, features = core_available, ncol = 4, pt.size = 0.3,
            order = TRUE) + plot_annotation(title = "Mouse: Core Cell Type Markers")
dev.off()

# --- A6. 基于 marker 手动注释 ---
cat("\n>>> 手动分配细胞类型...\n")

# Default: 使用 SingleR 注释
seu_mouse$cell_type <- seu_mouse$singleR_main

# 检查每个 cluster 中 marker 表达，精细调整
# (这里给出基于自动注释 + cluster 表达谱的逻辑框架)
# 用户可根据 DotPlot 结果进一步微调

# --- A7. 注释结果 UMAP ---
p_auto <- DimPlot(seu_mouse, group.by = "singleR_main", label = TRUE,
                  repel = TRUE, pt.size = 0.2) + ggtitle("Mouse: SingleR Annotation")
p_ct   <- DimPlot(seu_mouse, group.by = "cell_type", label = TRUE,
                  repel = TRUE, pt.size = 0.2) + ggtitle("Mouse: Cell Type")

pdf(here("GSE206861/results/mouse_annotation_umap.pdf"), width = 20, height = 8)
p_auto + p_ct
dev.off()

# --- A8. 细胞比例 ---
cell_prop <- seu_mouse@meta.data %>%
  group_by(orig.ident, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(proportion = n / sum(n) * 100)

pdf(here("GSE206861/results/mouse_cell_proportion.pdf"), width = 10, height = 6)
ggplot(cell_prop, aes(x = orig.ident, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Sample", y = "Proportion", fill = "Cell Type",
       title = "Mouse: Cell Type Proportion (CF vs Sham)") +
  scale_fill_viridis_d(option = "turbo") +
  theme_minimal()
dev.off()

cat("\n  小鼠注释完成。\n")

# ============================================================================
# Part B: 人类 ACP 样本注释
# ============================================================================

cat("\n========== 人类 ACP 注释 ==========\n")

# --- B1. 人类参考 ---
cat("\n>>> 加载人类参考...\n")
ref_human <- BlueprintEncodeData()

DefaultAssay(seu_acp) <- "SCT"
human_counts <- GetAssayData(seu_acp, layer = "counts")

pred_human <- SingleR(
  test   = human_counts,
  ref    = ref_human,
  labels = ref_human$label.main,
  de.method = "classic"
)

seu_acp$singleR_main <- pred_human$labels[match(colnames(seu_acp), rownames(pred_human))]

cat("  SingleR 注释结果:\n")
print(table(seu_acp$singleR_main))

# --- B2. 人类 marker ---
human_markers <- list(
  "Neuron"         = c("RBFOX3", "SYT1", "TUBB3", "SNAP25"),
  "Microglia"      = c("CX3CR1", "TMEM119", "P2RY12", "AIF1"),
  "Astrocyte"      = c("GFAP", "AQP4", "ALDH1L1", "SLC1A3"),
  "Oligo"          = c("MBP", "MOG", "OLIG1", "OLIG2"),
  "Endothelial"    = c("CLDN5", "PECAM1", "CDH5"),
  "Pericyte"       = c("PDGFRB", "RGS5", "ANPEP"),
  "T_cell"         = c("CD3D", "CD3E", "TRBC2"),
  "B_cell"         = c("CD79A", "MS4A1"),
  "Macrophage"     = c("CD68", "CD163", "CSF1R"),
  "Epithelial"     = c("EPCAM", "KRT8", "KRT18", "KRT19")  # ACP 是上皮来源肿瘤
)

# DotPlot
pdf(here("GSE206861/results/human_markers_dotplot.pdf"), width = 20, height = 8)
DotPlot(seu_acp, features = unique(unlist(human_markers)), assay = "SCT",
        group.by = "seurat_clusters") + RotatedAxis() +
  ggtitle("Human ACP: Marker Gene Expression by Cluster")
dev.off()

# --- B3. Feature Plot ---
human_core <- c("RBFOX3", "CX3CR1", "TMEM119", "GFAP", "MBP", "CD68",
                "EPCAM", "PECAM1", "CD3D", "MS4A1")
human_available <- intersect(human_core, rownames(seu_acp))

pdf(here("GSE206861/results/human_feature_core.pdf"), width = 16, height = 14)
FeaturePlot(seu_acp, features = human_available, ncol = 4, pt.size = 0.5,
            order = TRUE) + plot_annotation(title = "Human ACP: Core Cell Type Markers")
dev.off()

# --- B4. UMAP ---
pdf(here("GSE206861/results/human_annotation_umap.pdf"), width = 10, height = 8)
DimPlot(seu_acp, group.by = "singleR_main", label = TRUE, repel = TRUE,
        pt.size = 0.5) + ggtitle("Human ACP: SingleR Annotation")
dev.off()

# ============================================================================
# 保存
# ============================================================================

cat("\n========== 保存注释结果 ==========\n")
saveRDS(seu_mouse, here("GSE206861/data/processed/seu_mouse_annotated.rds"))
saveRDS(seu_acp,   here("GSE206861/data/processed/seu_acp_annotated.rds"))

cat("\n✓ 细胞注释完成！\n")
cat("  保存: data/processed/seu_mouse_annotated.rds\n")
cat("  保存: data/processed/seu_acp_annotated.rds\n")
cat("  接下来运行: 05_differential_expression.R\n")
