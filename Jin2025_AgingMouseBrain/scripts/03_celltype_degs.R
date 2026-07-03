# =============================================================================
# Phase 3: 分性别/分细胞类型 DEG分析
# =============================================================================
# anatomical set: Female Aged vs Adult, Male Aged vs Adult
# dissection set: Aged vs Adult (overall, due to sex imbalance)
# 对关键神经元亚类(ARH/DMH/PVH/SCH等)同样分析
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(EnhancedVolcano)
  library(ggrepel)
})

cat("\n============================================\n")
cat(" Phase 3: 差异表达分析\n")
cat("============================================\n\n")

DATA_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus"
RESULT_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/results"
dir.create(file.path(RESULT_DIR, "DEGs"), recursive = TRUE, showWarnings = FALSE)

# ---- 配置 ----
# 需要至少这么多细胞才做DEG
MIN_CELLS_PER_GROUP <- 10
# DEG参数
DEG_MIN_PCT <- 0.1
DEG_LOGFC <- 0.1

# ===========================================================================
# DEG分析函数
# ===========================================================================
run_degs <- function(obj, cell_types, ident_col, ident1, ident2,
                     subgroup_col = NULL, subgroup_val = NULL,
                     label_prefix = "") {
  # 对指定细胞类型运行 DEG 分析。
  # ident_col 定义比较组；subgroup_col/subgroup_val 可用于按性别或其他
  # 元数据字段进一步筛选细胞。
  results <- list()

  for (ct in cell_types) {
    cat(sprintf("\n--- %s ---\n", ct))

    # 子集
    cells_ct <- WhichCells(obj, expression = cell_type_major == ct)
    if (subgroup_col %in% colnames(obj@meta.data) && !is.null(subgroup_val)) {
      cells_sg <- WhichCells(obj, expression = !!sym(subgroup_col) == subgroup_val)
      cells_use <- intersect(cells_ct, cells_sg)
      sub_label <- sprintf("%s_%s", ct, subgroup_val)
    } else {
      cells_use <- cells_ct
      sub_label <- ct
    }

    if (length(cells_use) < MIN_CELLS_PER_GROUP * 2) {
      cat(sprintf("  ⚠️ 细胞数不足 (%d), 跳过\n", length(cells_use)))
      next
    }

    # 检查两组是否都有足够细胞
    sub_obj <- subset(obj, cells = cells_use)
    Idents(sub_obj) <- ident_col
    idents_present <- unique(Idents(sub_obj))

    if (!(ident1 %in% idents_present) || !(ident2 %in% idents_present)) {
      cat(sprintf("  ⚠️ %s或%s不存在, 跳过\n", ident1, ident2))
      next
    }

    n1 <- sum(Idents(sub_obj) == ident1)
    n2 <- sum(Idents(sub_obj) == ident2)
    cat(sprintf("  %s=%d, %s=%d\n", ident1, n1, ident2, n2))

    if (n1 < MIN_CELLS_PER_GROUP || n2 < MIN_CELLS_PER_GROUP) {
      cat(sprintf("  ⚠️ 组内细胞数不足, 跳过\n"))
      next
    }

    # 运行FindMarkers
    tryCatch({
      deg <- FindMarkers(
        sub_obj,
        ident.1 = ident1,
        ident.2 = ident2,
        min.pct = DEG_MIN_PCT,
        logfc.threshold = DEG_LOGFC,
        test.use = "wilcox"
      )

      deg$gene <- rownames(deg)
      deg <- deg[order(deg$avg_log2FC, decreasing = TRUE), ]

      n_sig <- sum(deg$p_val_adj < 0.05, na.rm = TRUE)
      n_up <- sum(deg$p_val_adj < 0.05 & deg$avg_log2FC > 0, na.rm = TRUE)
      n_down <- sum(deg$p_val_adj < 0.05 & deg$avg_log2FC < 0, na.rm = TRUE)

      cat(sprintf("  DEGs: %d sig (↑%d ↓%d)\n", n_sig, n_up, n_down))

      # 保存
      fname <- file.path(RESULT_DIR, "DEGs",
                         sprintf("%s%s_DEG.csv", label_prefix, gsub("/", "_", sub_label)))
      fwrite(deg, fname)

      # 火山图 (如果有显著DEG)
      if (n_sig >= 5) {
        deg$signif <- "Not Sig"
        deg$signif[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0.5] <- "Up"
        deg$signif[deg$p_val_adj < 0.05 & deg$avg_log2FC < -0.5] <- "Down"

        top_genes <- deg %>%
          filter(p_val_adj < 0.05) %>%
          arrange(desc(abs(avg_log2FC))) %>%
          head(15) %>%
          pull(gene)

        p <- ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300),
                             color = signif)) +
          geom_point(size = 0.5, alpha = 0.6) +
          scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5",
                                         "Not Sig" = "grey80")) +
          geom_text_repel(
            data = subset(deg, gene %in% top_genes),
            aes(label = gene), size = 3, max.overlaps = 30
          ) +
          geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.3) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
          ggtitle(sprintf("%s: %s vs %s", sub_label, ident1, ident2)) +
          theme_minimal(base_size = 11)

        ggsave(file.path(RESULT_DIR, "DEGs",
                         sprintf("%s%s_volcano.png", label_prefix, gsub("/", "_", sub_label))),
               p, width = 8, height = 7, dpi = 150)
      }

      results[[sub_label]] <- deg
      gc()

    }, error = function(e) {
      cat(sprintf("  ❌ DEG失败: %s\n", e$message))
    })
  }

  return(results)
}

# ===========================================================================
# 加载数据
# ===========================================================================
cat("\n[1] 加载Seurat对象...\n")
rds_anat <- file.path(DATA_DIR, "hypo_anatomical_seurat.rds")
rds_diss <- file.path(DATA_DIR, "hypo_dissection_seurat.rds")

# 优先使用anatomical (性别均衡)
if (file.exists(rds_anat)) {
  hypo_anat <- readRDS(rds_anat)
  cat(sprintf("  anatomical: %d cells\n", ncol(hypo_anat)))
} else {
  cat("  ⚠️ anatomical RDS不存在, 从h5ad加载...\n")
}

if (file.exists(rds_diss)) {
  hypo_diss <- readRDS(rds_diss)
  cat(sprintf("  dissection: %d cells\n", ncol(hypo_diss)))
} else {
  cat("  ⚠️ dissection RDS不存在\n")
}

# ===========================================================================
# 获取cell_type列表
# ===========================================================================
all_ct <- unique(hypo_anat$cell_type_major)
# 过滤掉细胞数太少的类型
ct_counts <- table(hypo_anat$cell_type_major)
major_ct <- names(ct_counts[ct_counts >= 20])
minor_ct <- names(ct_counts[ct_counts < 20])

cat(sprintf("\n  主要细胞类型(>=20 cells): %d\n", length(major_ct)))
cat(sprintf("  次要细胞类型(<20 cells): %d (跳过)\n", length(minor_ct)))

# ===========================================================================
# 2. Anatomical set: 分性别DEG
# ===========================================================================
cat("\n[2] Anatomical set: Female Aged vs Adult...\n")
results_anat_F <- run_degs(
  hypo_anat, major_ct,
  ident_col = "age_group", ident1 = "Aged(18m)", ident2 = "Adult(2m)",
  subgroup_col = "sex_group", subgroup_val = "Female",
  label_prefix = "anat_F_"
)

cat("\n[3] Anatomical set: Male Aged vs Adult...\n")
results_anat_M <- run_degs(
  hypo_anat, major_ct,
  ident_col = "age_group", ident1 = "Aged(18m)", ident2 = "Adult(2m)",
  subgroup_col = "sex_group", subgroup_val = "Male",
  label_prefix = "anat_M_"
)

cat("\n[4] Anatomical set: Overall Aged vs Adult (不分子组)...\n")
results_anat_all <- run_degs(
  hypo_anat, major_ct,
  ident_col = "age_group", ident1 = "Aged(18m)", ident2 = "Adult(2m)",
  label_prefix = "anat_all_"
)

# ===========================================================================
# 3. 神经元亚类 DEG (如果neuronal_subclass存在)
# ===========================================================================
if ("neuronal_subclass" %in% colnames(hypo_anat@meta.data)) {
  cat("\n[5] 关键神经元亚类 DEG...\n")

  key_subclasses <- c("ARH_Neuron", "DMH_Neuron", "PVH_Neuron",
                      "VMH_Neuron", "SCH_Neuron", "LHA_Neuron",
                      "MPO_Neuron", "TU_Neuron", "MM_Neuron")

  # 先查看实际有哪些亚类
  actual_subclasses <- unique(hypo_anat$neuronal_subclass)
  cat(sprintf("  实际神经元亚类: %s\n", paste(actual_subclasses, collapse=", ")))

  use_subclasses <- intersect(key_subclasses, actual_subclasses)
  cat(sprintf("  将分析 %d 个亚类\n", length(use_subclasses)))

  for (sc in use_subclasses) {
    cat(sprintf("\n  --- %s ---\n", sc))
    cells_sc <- WhichCells(hypo_anat, expression = neuronal_subclass == sc)
    cat(sprintf("  Total cells: %d\n", length(cells_sc)))

    if (length(cells_sc) < 20) {
      cat("  ⚠️ 细胞数不足, 跳过\n")
      next
    }

    sub_obj <- subset(hypo_anat, cells = cells_sc)

    # 分性别
    for (sex in c("Female", "Male")) {
      cells_sex <- WhichCells(sub_obj, expression = sex_group == sex)
      if (length(cells_sex) < 20) next

      sub_sex <- subset(sub_obj, cells = cells_sex)
      Idents(sub_sex) <- "age_group"

      if (sum(Idents(sub_sex) == "Aged(18m)") < MIN_CELLS_PER_GROUP ||
          sum(Idents(sub_sex) == "Adult(2m)") < MIN_CELLS_PER_GROUP) next

      tryCatch({
        deg <- FindMarkers(sub_sex, ident.1 = "Aged(18m)", ident.2 = "Adult(2m)",
                           min.pct = DEG_MIN_PCT, logfc.threshold = DEG_LOGFC,
                           test.use = "wilcox")
        deg$gene <- rownames(deg)
        deg <- deg[order(deg$avg_log2FC, decreasing = TRUE), ]

        n_sig <- sum(deg$p_val_adj < 0.05, na.rm = TRUE)
        cat(sprintf("    %s %s: %d DEGs\n", sc, sex, n_sig))

        fname <- file.path(RESULT_DIR, "DEGs",
                           sprintf("subclass_%s_%s_DEG.csv", sc, sex))
        fwrite(deg, fname)
      }, error = function(e) {
        cat(sprintf("    ❌ %s\n", e$message))
      })
    }
  }
}

# ===========================================================================
# 4. Dissection set: 仅总体DEG (性别失衡)
# ===========================================================================
if (exists("hypo_diss")) {
  cat("\n[6] Dissection set: Overall Aged vs Adult...\n")

  diss_ct <- names(which(table(hypo_diss$cell_type_major) >= 20))

  results_diss <- run_degs(
    hypo_diss, diss_ct,
    ident_col = "age_group", ident1 = "Aged(18m)", ident2 = "Adult(2m)",
    label_prefix = "diss_all_"
  )
}

# ===========================================================================
# 5. DEG数量汇总
# ===========================================================================
cat("\n[7] 生成DEG汇总表...\n")

collect_degs <- function(result_list, label) {
  summary <- data.frame(
    cell_type = character(),
    n_sig = integer(),
    n_up = integer(),
    n_down = integer(),
    stringsAsFactors = FALSE
  )
  for (ct in names(result_list)) {
    deg <- result_list[[ct]]
    if (is.null(deg)) next
    summary <- rbind(summary, data.frame(
      cell_type = ct,
      n_sig = sum(deg$p_val_adj < 0.05, na.rm = TRUE),
      n_up  = sum(deg$p_val_adj < 0.05 & deg$avg_log2FC > 0, na.rm = TRUE),
      n_down = sum(deg$p_val_adj < 0.05 & deg$avg_log2FC < 0, na.rm = TRUE)
    ))
  }
  summary <- summary[order(summary$n_sig, decreasing = TRUE), ]
  fwrite(summary, file.path(RESULT_DIR, "DEGs", sprintf("DEG_summary_%s.csv", label)))
  return(summary)
}

collect_degs(results_anat_F, "anatomical_Female")
collect_degs(results_anat_M, "anatomical_Male")
collect_degs(results_anat_all, "anatomical_overall")

if (exists("results_diss")) {
  collect_degs(results_diss, "dissection_overall")
}

cat("\n============================================\n")
cat(" Phase 3 完成!\n")
cat("============================================\n")
cat(sprintf("  结果目录: %s/DEGs/\n", RESULT_DIR))
