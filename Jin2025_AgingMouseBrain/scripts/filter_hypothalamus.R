# =============================================================================
# 从全脑数据中提取下丘脑 (Hypothalamus) 区域数据
# =============================================================================
# 基于 Allen CCFv3 脑区注释，筛选下丘脑及相关亚区
#
# 下丘脑在CCFv3中的结构ID:
#   HY (Hypothalamus): 1097
#   主要子区域:
#     ARH - Arcuate hypothalamic nucleus
#     DMH - Dorsomedial nucleus of the hypothalamus
#     LHA - Lateral hypothalamic area
#     MPN - Medial preoptic nucleus
#     PH  - Posterior hypothalamic nucleus
#     PVH - Paraventricular hypothalamic nucleus
#     PVp - Periventricular hypothalamic nucleus, posterior part
#     SCH - Suprachiasmatic nucleus
#     SO  - Supraoptic nucleus
#     VMH - Ventromedial hypothalamic nucleus
#     TU  - Tuberal nucleus
#     ...
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# 配置 ----
PROJECT_DIR <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
if (PROJECT_DIR == "" || is.null(PROJECT_DIR)) {
  PROJECT_DIR <- getwd()
}
DATA_DIR <- file.path(PROJECT_DIR, "data")
OUT_DIR <- file.path(DATA_DIR, "hypothalamus")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat(" 提取下丘脑区域数据\n")
cat("============================================\n\n")

# ============================================================================
# 1. 从单细胞数据中筛选下丘脑细胞 ----
# ============================================================================
cat("[1/3] 从 scRNA-seq 数据中提取下丘脑细胞...\n")

cell_meta_file <- file.path(DATA_DIR, "single_cell", "metadata", "cell_metadata.csv")
cluster_annot_file <- file.path(DATA_DIR, "single_cell", "metadata",
                                "cell_cluster_annotations.csv")
mapping_file <- file.path(DATA_DIR, "single_cell", "taxonomy",
                          "cell_cluster_mapping_annotations.csv")

if (!file.exists(cell_meta_file)) {
  cat(sprintf("  ✗ cell_metadata.csv 未找到! 请先运行 01_download_metadata.sh\n"))
  quit(status = 1)
}

cat("  读取 cell_metadata.csv (文件较大，请耐心等待)...\n")
cell_meta <- fread(cell_meta_file, showProgress = TRUE)
cat(sprintf("  总细胞数: %s\n", format(nrow(cell_meta), big.mark = ",")))

# 检查哪些列可用于筛选脑区
cat("\n  可用列名:\n")
cat(paste("   ", colnames(cell_meta), collapse = "\n    "), "\n")

# 查找与脑区/解剖结构相关的列
region_cols <- grep("parcellation|region|structure|anatomy|area|location",
                    colnames(cell_meta), value = TRUE, ignore.case = TRUE)
cat(sprintf("\n  可能的脑区相关列: %s\n",
            paste(region_cols, collapse = ", ")))

# 查找下丘脑细胞
# 策略1: 在parcellation/region列中搜索 "hypothalamus" / "HY" 关键词
# 策略2: 使用 CCF 结构ID (1097 = hypothalamus)
# 策略3: 使用 subclass/class/tissue 等聚类注释

hypo_cells <- NULL
hypo_method <- ""

# 尝试策略1: 文本匹配
text_cols <- c("parcellation_label", "parcellation_name", "region_label",
               "structure_label", "anatomy", "tissue", "brain_region")

for (col in text_cols) {
  if (col %in% colnames(cell_meta)) {
    hypo_idx <- grepl("hypothalamus|hypothal|HY[^P]|Arcuate|ARH|DMH|LHA|MPN|PVH|VMH|SCH|SO[^M]|Tuberal",
                      cell_meta[[col]], ignore.case = TRUE)
    if (sum(hypo_idx) > 0) {
      cat(sprintf("\n  在列 '%s' 中找到 %s 个下丘脑相关细胞\n",
                  col, format(sum(hypo_idx), big.mark = ",")))
      if (is.null(hypo_cells)) {
        hypo_cells <- which(hypo_idx)
        hypo_method <- col
      } else {
        hypo_cells <- union(hypo_cells, which(hypo_idx))
        hypo_method <- paste(hypo_method, col, sep = "+")
      }
    }
  }
}

# 尝试策略2: CCF结构ID匹配
id_cols <- c("parcellation_id", "structure_id", "region_id", "ccf_id")

for (col in id_cols) {
  if (col %in% colnames(cell_meta)) {
    # CCFv3: 1097 = Hypothalamus
    # 子区域ID范围族
    hypo_ids <- cell_meta[[col]] == 1097 |
      cell_meta[[col]] %in% c(185, 223, 249, 272, 315, 318, 321, 324, 327,
                              330, 349, 380, 400, 560, 613, 623, 629, 635,
                              648, 673, 684, 700, 704, 723, 732, 746)
    if (sum(hypo_ids) > 0) {
      cat(sprintf("\n  在列 '%s' 中通过CCF ID找到 %s 个下丘脑细胞\n",
                  col, format(sum(hypo_ids), big.mark = ",")))
      if (is.null(hypo_cells)) {
        hypo_cells <- which(hypo_ids)
        hypo_method <- paste0(col, "(CCF ID)")
      } else {
        hypo_cells <- union(hypo_cells, which(hypo_ids))
      }
    }
  }
}

# 尝试策略3: 从cluster/子类注释中匹配
if (!is.null(hypo_cells)) {
  cat(sprintf("\n  ✅ 共筛选出 %s 个下丘脑细胞\n",
              format(length(hypo_cells), big.mark = ",")))
  cat(sprintf("     方法: %s\n", hypo_method))

  hypo_meta <- cell_meta[hypo_cells, ]

  # 保存筛选后的metadata
  hypo_meta_file <- file.path(OUT_DIR, "hypothalamus_cell_metadata.csv")
  fwrite(hypo_meta, hypo_meta_file)
  cat(sprintf("     已保存: %s\n", hypo_meta_file))
} else {
  cat("\n  ⚠️ 未找到直接的下丘脑标记列!\n")
  cat("     将保存完整cell_metadata供手动探索\n")
  cat("     请在ABC Atlas Explorer中确认脑区注释列名\n")

  # 查看数据中的parcellation值分布
  for (col in grep("parcellation|region|area|tissue|structure",
                   colnames(cell_meta), value = TRUE, ignore.case = TRUE)) {
    cat(sprintf("\n  列 '%s' 的唯一值 (前30个):\n", col))
    vals <- unique(cell_meta[[col]])
    print(head(vals, 30))
  }
}

# ============================================================================
# 2. 从 MERFISH 空间数据中筛选下丘脑坐标 ----
# ============================================================================
cat("\n[2/3] 从 MERFISH 空间数据中提取下丘脑坐标...\n")

ccf_file <- file.path(DATA_DIR, "spatial", "MERFISH_638850", "metadata",
                      "ccf_coordinates.csv")

if (file.exists(ccf_file)) {
  cat("  读取 ccf_coordinates.csv...\n")
  ccf <- fread(ccf_file, showProgress = TRUE)
  cat(sprintf("  总细胞数 (MERFISH): %s\n", format(nrow(ccf), big.mark = ",")))
  cat(sprintf("  可用列: %s\n", paste(colnames(ccf), collapse = ", ")))

  # MERFISH ccf_coordinates 通常包含:
  # cell_label, x, y, z, parcellation_label, parcellation_id 等

  # 搜索下丘脑
  if ("parcellation_label" %in% colnames(ccf)) {
    hypo_ccf_idx <- grepl("hypothalamus|hypothal|HY[^P]",
                          ccf$parcellation_label, ignore.case = TRUE)
  } else if ("parcellation_name" %in% colnames(ccf)) {
    hypo_ccf_idx <- grepl("hypothalamus|hypothal|HY[^P]",
                          ccf$parcellation_name, ignore.case = TRUE)
  } else if ("parcellation_id" %in% colnames(ccf)) {
    # 需要CCF结构树来获取下丘脑的子ID
    cat("  ⚠️ 仅有parcellation_id，需要CCF结构树来解析。\n")
    cat("     下载Allen CCF结构树进行匹配...\n")

    # 尝试下载Allen CCF structure tree
    structure_tree_url <- "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/WMB-CCF-structure-tree/20230630/structure_tree.csv"
    structure_tree_file <- file.path(DATA_DIR, "spatial", "structure_tree.csv")

    tryCatch({
      if (!file.exists(structure_tree_file)) {
        download.file(structure_tree_url, structure_tree_file, mode = "wb")
      }
      structure_tree <- fread(structure_tree_file)
      # 找到下丘脑及其所有子结构
      # CCFv3中下丘脑 (HY) 的ID是1097
      hypo_structure_ids <- subset(structure_tree,
                                   grepl("hypothalamus|hypothal",
                                         name, ignore.case = TRUE) |
                                   id == 1097 |
                                   structure_set_ids == 1097)
      cat(sprintf("     下丘脑相关CCF结构数: %d\n", nrow(hypo_structure_ids)))

      # 还需要递归获取所有子结构...
      hypo_ccf_idx <- ccf$parcellation_id %in% c(1097, hypo_structure_ids$id)
    }, error = function(e) {
      cat(sprintf("     下载结构树失败: %s\n", e$message))
      hypo_ccf_idx <- ccf$parcellation_id == 1097
    })
  } else {
    hypo_ccf_idx <- rep(FALSE, nrow(ccf))
    cat("  ⚠️ 未找到parcellation相关列\n")
  }

  if (sum(hypo_ccf_idx) > 0) {
    cat(sprintf("  ✅ MERFISH下丘脑细胞数: %s\n",
                format(sum(hypo_ccf_idx), big.mark = ",")))

    hypo_ccf <- ccf[hypo_ccf_idx, ]
    hypo_ccf_file <- file.path(OUT_DIR, "hypothalamus_MERFISH_coordinates.csv")
    fwrite(hypo_ccf, hypo_ccf_file)
    cat(sprintf("     已保存: %s\n", hypo_ccf_file))

    # 保存下丘脑细胞的cell_label，用于后续从表达矩阵中提取
    hypo_cell_labels <- hypo_ccf$cell_label
    cell_label_file <- file.path(OUT_DIR, "hypothalamus_MERFISH_cell_labels.csv")
    fwrite(data.table(cell_label = hypo_cell_labels), cell_label_file)
    cat(sprintf("     已保存细胞标签: %s\n", cell_label_file))
  } else {
    cat("  ⚠️ 未找到明确的下丘脑标记\n")
  }
} else {
  cat(sprintf("  ✗ ccf_coordinates.csv 未找到: %s\n", ccf_file))
  cat("     请先运行 01_download_metadata.sh\n")
}

# ============================================================================
# 3. 提取下丘脑年龄差异基因 ----
# ============================================================================
cat("\n[3/3] 提取下丘脑相关年龄差异基因...\n")

deg_file <- file.path(DATA_DIR, "single_cell", "taxonomy", "aging_degenes.csv")

if (file.exists(deg_file)) {
  deg <- fread(deg_file)
  cat(sprintf("  年龄差异基因总数: %s\n", format(nrow(deg), big.mark = ",")))
  cat(sprintf("  可用列: %s\n", paste(colnames(deg), collapse = ", ")))

  # 如果deg文件包含cluster/subclass信息，筛选与下丘脑相关的
  # 如果deg文件包含区域信息，直接筛选
  hypo_deg <- deg

  region_in_deg <- grep("region|parcellation|area|structure|cluster|subclass|class|tissue",
                        colnames(deg), value = TRUE, ignore.case = TRUE)

  if (length(region_in_deg) > 0) {
    cat(sprintf("\n  区域相关列: %s\n", paste(region_in_deg, collapse = ", ")))
    # 检查是否有与下丘脑聚类的关联
    for (col in region_in_deg) {
      hypo_match <- grepl("hypothalamus|hypothal|tanycyte|ependymal|ARC|DMH|PVH|VMH",
                          deg[[col]], ignore.case = TRUE)
      if (sum(hypo_match) > 0) {
        cat(sprintf("  在 '%s' 中发现 %d 个下丘脑相关记录\n",
                    col, sum(hypo_match)))
      }
    }
  }

  hypo_deg_file <- file.path(OUT_DIR, "hypothalamus_related_aging_degenes.csv")
  fwrite(hypo_deg, hypo_deg_file)
  cat(sprintf("  已保存年龄差异基因: %s\n", hypo_deg_file))
} else {
  cat(sprintf("  ✗ aging_degenes.csv 未找到: %s\n", deg_file))
}

# ============================================================================
# 4. 总结 ----
# ============================================================================
cat("\n============================================\n")
cat(" 下丘脑数据提取完成!\n")
cat("============================================\n")
cat(sprintf(" 输出目录: %s\n", OUT_DIR))
cat("\n")
cat(" 生成的文件:\n")
for (f in list.files(OUT_DIR, full.names = TRUE)) {
  cat(sprintf("   %s (%s)\n", basename(f),
              format(file.info(f)$size, big.mark = ",")))
}
cat("\n")
cat(" 注意事项:\n")
cat("  1. 如果自动筛选不理想，请根据cell_metadata中的实际列名手动调整\n")
cat("  2. 可访问ABC Atlas Explorer验证脑区注释:\n")
cat("     https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas\n")
cat("  3. CCFv3完整结构树见:\n")
cat("     https://alleninstitute.github.io/abc_atlas_access/descriptions/CCF_parcellation_annotations.html\n")
