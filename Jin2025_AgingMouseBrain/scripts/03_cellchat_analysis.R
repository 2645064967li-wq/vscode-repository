# =============================================================================
# CellChat Analysis: Young vs Aged 小鼠下丘脑细胞通讯比较
# Jin et al. 2025 Nature
# 聚焦: Microglia ↔ 其他细胞类群的配体-受体互作
# =============================================================================
set.seed(42)

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB for parallel

cat("\n============================================\n")
cat(" CellChat 衰老小鼠下丘脑细胞通讯分析\n")
cat("============================================\n\n")

# ---- 0. Load packages ----
suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
})
cat("Packages loaded\n")

# ---- 1. Paths ----
DATA_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/data/cellchat"
RESULT_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/results/cellchat"
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULT_DIR, "microglia"), recursive = TRUE, showWarnings = FALSE)

# ---- 2. Load data ----
cat("\n[Step 1] Loading data...\n")

young_meta <- fread(file.path(DATA_DIR, "young_metadata.csv"))
aged_meta  <- fread(file.path(DATA_DIR, "aged_metadata.csv"))

# Filter out "Unclassified" and tiny cell types
MIN_CELLS <- 20
young_meta <- young_meta[young_meta$cell_type != "Unclassified", ]
aged_meta  <- aged_meta[aged_meta$cell_type != "Unclassified", ]

ct_counts_y <- young_meta[, .N, by = cell_type][N >= MIN_CELLS]
ct_counts_a <- aged_meta[, .N, by = cell_type][N >= MIN_CELLS]
keep_ct <- intersect(ct_counts_y$cell_type, ct_counts_a$cell_type)

young_meta <- young_meta[young_meta$cell_type %in% keep_ct, ]
aged_meta  <- aged_meta[aged_meta$cell_type %in% keep_ct, ]

cat(sprintf("Young: %d cells, %d cell types\n", nrow(young_meta), length(keep_ct)))
cat(sprintf("Aged:  %d cells, %d cell types\n", nrow(aged_meta), length(keep_ct)))
cat("Cell types:", paste(sort(keep_ct), collapse = ", "), "\n")

# ---- 3. Create CellChat objects ----
cat("\n[Step 2] Creating CellChat objects...\n")

for (label in c("young", "aged")) {
  cat(sprintf("\n  --- %s ---\n", label))

  counts_file <- file.path(DATA_DIR, paste0(label, "_raw_counts.csv.gz"))
  meta <- if (label == "young") young_meta else aged_meta

  # Load counts (genes x cells)
  cat("    Loading counts (may take a few minutes)...\n")
  counts_dt <- fread(counts_file)
  gene_names <- counts_dt[[1]]
  counts_mat <- as.matrix(counts_dt[, -1, drop = FALSE])
  rownames(counts_mat) <- gene_names

  # Filter to cells in metadata
  valid_cells <- intersect(colnames(counts_mat), meta$cell_label)
  counts_mat <- counts_mat[, valid_cells, drop = FALSE]
  meta_use <- meta[match(valid_cells, meta$cell_label), ]

  cat(sprintf("    Matrix: %d genes x %d cells\n", nrow(counts_mat), ncol(counts_mat)))

  # Create CellChat
  cc <- createCellChat(object = counts_mat, meta = meta_use, group.by = "cell_type")
  cc@DB <- CellChatDB.mouse
  cc <- subsetData(cc)

  # subsetData reduces to signaling genes only
  cat(sprintf("    Signaling genes: %d\n", nrow(cc@data.signaling)))

  # Store
  assign(paste0("cc_", label), cc)
  rm(counts_dt, counts_mat)
  gc()
}

# ---- 4. Run CellChat pipeline ----
cat("\n[Step 3] Running CellChat pipeline...\n")

run_chat <- function(cc, label) {
  cat(sprintf("\n  === %s ===\n", label))

  cat("  identifyOverExpressedGenes...\n")
  cc <- identifyOverExpressedGenes(cc)
  cat("  identifyOverExpressedInteractions...\n")
  cc <- identifyOverExpressedInteractions(cc)
  cat("  computeCommunProb (truncatedMean)...\n")
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                           raw.use = TRUE, population.size = TRUE)
  cat("  filterCommunication (min.cells=10)...\n")
  cc <- filterCommunication(cc, min.cells = 10)

  cat("  computeCommunProbPathway...\n")
  cc <- tryCatch({
    computeCommunProbPathway(cc)
  }, error = function(e) {
    cat("    computeCommunProbPathway failed, building manually:", conditionMessage(e), "\n")
    # Manually construct netP from net$prob
    prob_mat <- cc@net$prob
    ct_names <- levels(cc@idents)
    n_ct <- length(ct_names)
    cc@netP$pathways <- "All"
    cc@netP$prob <- array(prob_mat, dim = c(n_ct, n_ct, 1),
                          dimnames = list(ct_names, ct_names, "All"))
    cc
  })

  cat("  aggregateNet...\n")
  cc <- aggregateNet(cc)

  return(cc)
}

cc_young <- run_chat(cc_young, "Young(2m)")
cc_aged  <- run_chat(cc_aged, "Aged(18m)")

# Save
cat("\n  Saving CellChat objects...\n")
saveRDS(cc_young, file.path(RESULT_DIR, "cellchat_young.rds"))
saveRDS(cc_aged,  file.path(RESULT_DIR, "cellchat_aged.rds"))

# ---- 5. Merge and compare ----
cat("\n[Step 4] Merging and comparing...\n")

common_ct <- intersect(levels(cc_young@idents), levels(cc_aged@idents))
cat(sprintf("  Common cell types: %d\n", length(common_ct)))

# Sync netP structures before merge
cat("  Syncing netP structures...\n")
ypaths <- cc_young@netP$pathways
apaths <- cc_aged@netP$pathways
cat(sprintf("    Young pathways: %d, Aged pathways: %d\n", length(ypaths), length(apaths)))

if (length(ypaths) > 0 && length(apaths) == 1 && apaths[1] == "All") {
  cat("    Aged has fallback netP. Adapting to match Young structure...\n")
  ct_names <- levels(cc_aged@idents)
  n_ct <- length(ct_names)
  # Create a prob array with same number of pathways as Young, all zero
  n_paths <- length(ypaths)
  cc_aged@netP$pathways <- ypaths
  cc_aged@netP$prob <- array(0, dim = c(n_ct, n_ct, n_paths),
                             dimnames = list(ct_names, ct_names, ypaths))
}


cc_young_sub <- subsetCellChat(cc_young, idents.use = common_ct)
cc_aged_sub  <- subsetCellChat(cc_aged, idents.use = common_ct)

# Check if both have valid netP
has_pathways <- length(cc_young_sub@netP$pathways) > 0 && length(cc_aged_sub@netP$pathways) > 0

if (has_pathways) {
  cc_merged <- mergeCellChat(list(cc_young_sub, cc_aged_sub),
                             cell.names = c("Young(2m)", "Aged(18m)"))
  saveRDS(cc_merged, file.path(RESULT_DIR, "cellchat_merged.rds"))
  cat("  Merged object saved\n")
} else {
  cat("  WARNING: Cannot merge. Saving individual objects only.\n")
  cc_merged <- NULL
  saveRDS(cc_young_sub, file.path(RESULT_DIR, "cellchat_young_sub.rds"))
  saveRDS(cc_aged_sub, file.path(RESULT_DIR, "cellchat_aged_sub.rds"))
}

# ---- 6. Overall comparison figures ----
cat("\n[Step 5] Overall comparison figures...\n")
fig_dir <- file.path(RESULT_DIR, "figures")

if (!is.null(cc_merged)) {
  # 6a. Interaction counts comparison
  pdf(file.path(fig_dir, "01_interaction_counts.pdf"), width = 8, height = 6)
  print(compareInteractions(cc_merged, show.legend = FALSE, group = c(1, 2)))
  dev.off()
  cat("  OK: 01_interaction_counts.pdf\n")

  # 6b. Differential interaction numbers
  pdf(file.path(fig_dir, "02_diff_interactions.pdf"), width = 8, height = 6)
  netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "count")
  dev.off()
  cat("  OK: 02_diff_interactions.pdf\n")

  # 6c. Differential interaction strength
  pdf(file.path(fig_dir, "03_diff_strength.pdf"), width = 8, height = 6)
  netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "weight")
  dev.off()
  cat("  OK: 03_diff_strength.pdf\n")

  # 6d. Signaling pathway rank comparison
  pdf(file.path(fig_dir, "04_pathway_rank.pdf"), width = 10, height = 12)
  rankNet(cc_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE,
          font.size = 8, title = "Pathway Information Flow: Young vs Aged")
  dev.off()
  cat("  OK: 04_pathway_rank.pdf\n")

  # 6e. Heatmap of overall signaling
  pdf(file.path(fig_dir, "05_signaling_heatmap.pdf"), width = 14, height = 7)
  netVisual_heatmap(cc_merged, comparison = c(1, 2), font.size = 7)
  dev.off()
  cat("  OK: 05_signaling_heatmap.pdf\n")

  # 6f. Bubble plot
  tryCatch({
    pdf(file.path(fig_dir, "06_signaling_bubble.pdf"), width = 16, height = 10)
    netVisual_bubble(cc_young_sub, cc_aged_sub, comparison = TRUE,
                     angle.x = 45, font.size = 8)
    dev.off()
    cat("  OK: 06_signaling_bubble.pdf\n")
  }, error = function(e) {
    cat("  Bubble plot skipped:", e$message, "\n")
  })
} else {
  cat("  Skipping merge-dependent figures (no merged object)\n")
}

# ---- 7. Microglia-focused analysis ----
cat("\n[Step 6] Microglia-focused analysis...\n")
mg_dir <- file.path(RESULT_DIR, "microglia")

if (!"Microglia" %in% common_ct) {
  cat("  WARNING: Microglia not found!\n")
} else {
  other_ct <- setdiff(common_ct, "Microglia")
  mg_idx_y <- which(levels(cc_young_sub@idents) == "Microglia")
  mg_idx_a <- which(levels(cc_aged_sub@idents) == "Microglia")
  other_idx_y <- which(levels(cc_young_sub@idents) %in% other_ct)
  other_idx_a <- which(levels(cc_aged_sub@idents) %in% other_ct)

  # 7a. Circle plots: Microglia outgoing
  cat("  7a. MG outgoing circle plots...\n")
  pdf(file.path(mg_dir, "MG01_outgoing_circle.pdf"), width = 14, height = 7)
  par(mfrow = c(1, 2))
  netVisual_circle(cc_young_sub@net$count,
                   sources.use = mg_idx_y,
                   targets.use = other_idx_y,
                   title.name = "Young(2m): MG → Others",
                   vertex.label.cex = 0.7, remove.isolate = TRUE)
  netVisual_circle(cc_aged_sub@net$count,
                   sources.use = mg_idx_a,
                   targets.use = other_idx_a,
                   title.name = "Aged(18m): MG → Others",
                   vertex.label.cex = 0.7, remove.isolate = TRUE)
  dev.off()
  cat("    OK\n")

  # 7b. Circle plots: Microglia incoming
  cat("  7b. MG incoming circle plots...\n")
  pdf(file.path(mg_dir, "MG02_incoming_circle.pdf"), width = 14, height = 7)
  par(mfrow = c(1, 2))
  netVisual_circle(cc_young_sub@net$count,
                   sources.use = other_idx_y,
                   targets.use = mg_idx_y,
                   title.name = "Young(2m): Others → MG",
                   vertex.label.cex = 0.7, remove.isolate = TRUE)
  netVisual_circle(cc_aged_sub@net$count,
                   sources.use = other_idx_a,
                   targets.use = mg_idx_a,
                   title.name = "Aged(18m): Others → MG",
                   vertex.label.cex = 0.7, remove.isolate = TRUE)
  dev.off()
  cat("    OK\n")

  # 7c. Signaling role scatter
  cat("  7c. Signaling role scatter...\n")
  pdf(file.path(mg_dir, "MG03_signaling_role.pdf"), width = 14, height = 7)
  par(mfrow = c(1, 2))
  netAnalysis_signalingRole_scatter(cc_young_sub,
    title = "Young(2m): Signaling Roles", font.size = 8)
  netAnalysis_signalingRole_scatter(cc_aged_sub,
    title = "Aged(18m): Signaling Roles", font.size = 8)
  dev.off()
  cat("    OK\n")

  # 7d. Microglia pathway changes
  if (!is.null(cc_merged)) {
    cat("  7d. MG pathway change ranking...\n")
    pathways <- cc_merged@netP$pathways
    mg_changes <- data.frame(pathway = pathways, stringsAsFactors = FALSE)

    for (i in seq_along(pathways)) {
      pw <- pathways[i]
      mg_changes$young_outgoing[i] <- if (pw %in% names(cc_young_sub@netP$prob))
        sum(cc_young_sub@netP$prob[[pw]][mg_idx_y, ]) else NA
      mg_changes$aged_outgoing[i] <- if (pw %in% names(cc_aged_sub@netP$prob))
        sum(cc_aged_sub@netP$prob[[pw]][mg_idx_a, ]) else NA
      mg_changes$young_incoming[i] <- if (pw %in% names(cc_young_sub@netP$prob))
        sum(cc_young_sub@netP$prob[[pw]][, mg_idx_y]) else NA
      mg_changes$aged_incoming[i] <- if (pw %in% names(cc_aged_sub@netP$prob))
        sum(cc_aged_sub@netP$prob[[pw]][, mg_idx_a]) else NA
    }

    mg_changes$outgoing_log2FC <- log2(
      (mg_changes$aged_outgoing + 1e-6) / (mg_changes$young_outgoing + 1e-6))
    mg_changes$incoming_log2FC <- log2(
      (mg_changes$aged_incoming + 1e-6) / (mg_changes$young_incoming + 1e-6))
    mg_changes <- mg_changes[order(abs(mg_changes$outgoing_log2FC), decreasing = TRUE), ]

    fwrite(mg_changes, file.path(mg_dir, "MG04_microglia_pathway_changes.csv"))

    cat("\n  Top pathways with MG outgoing INCREASE in aged:\n")
    top_up <- head(mg_changes[order(mg_changes$outgoing_log2FC, decreasing = TRUE),
                              c("pathway", "young_outgoing", "aged_outgoing", "outgoing_log2FC")], 10)
    print(top_up, row.names = FALSE)

    cat("\n  Top pathways with MG outgoing DECREASE in aged:\n")
    top_down <- head(mg_changes[order(mg_changes$outgoing_log2FC),
                                c("pathway", "young_outgoing", "aged_outgoing", "outgoing_log2FC")], 10)
    print(top_down, row.names = FALSE)
  } else {
    cat("  7d. Skipping MG pathway change ranking (no merged object)\n")
  }

  # 7e. Chord diagram for Microglia
  cat("  7e. MG chord diagrams...\n")
  tryCatch({
    pdf(file.path(mg_dir, "MG05_chord_outgoing.pdf"), width = 14, height = 7)
    par(mfrow = c(1, 2))
    netVisual_chord_gene(cc_young_sub, sources.use = mg_idx_y,
      targets.use = other_idx_y,
      title.name = "Young(2m): MG Outgoing",
      legend.pos.x = 8, legend.pos.y = 20)
    netVisual_chord_gene(cc_aged_sub, sources.use = mg_idx_a,
      targets.use = other_idx_a,
      title.name = "Aged(18m): MG Outgoing",
      legend.pos.x = 8, legend.pos.y = 20)
    dev.off()
    cat("    OK\n")
  }, error = function(e) {
    cat("    Chord diagram failed:", e$message, "\n")
  })
}

# ---- 8. Summary stats ----
cat("\n[Step 7] Saving summary statistics...\n")

# Per cell type network stats
for (label in c("young", "aged")) {
  cc <- if (label == "young") cc_young_sub else cc_aged_sub
  net_sum <- data.frame(
    cell_type = levels(cc@idents),
    outgoing_count = colSums(cc@net$count),
    incoming_count = rowSums(cc@net$count),
    outgoing_weight = colSums(cc@net$weight),
    incoming_weight = rowSums(cc@net$weight)
  )
  fwrite(net_sum, file.path(RESULT_DIR, paste0(label, "_net_stats.csv")))
}

# Summary
summary_df <- data.frame(
  Group = c("Young(2m)", "Aged(18m)"),
  Interactions = c(sum(cc_young_sub@net$count), sum(cc_aged_sub@net$count)),
  TotalWeight = c(sum(cc_young_sub@net$weight), sum(cc_aged_sub@net$weight)),
  Pathways = c(length(cc_young_sub@netP$pathways), length(cc_aged_sub@netP$pathways))
)
fwrite(summary_df, file.path(RESULT_DIR, "summary.csv"))
print(summary_df)

# ---- 9. Done ----
cat("\n============================================\n")
cat(" CellChat analysis complete!\n")
cat("============================================\n")
cat(sprintf(" Results: %s/\n", RESULT_DIR))
cat(sprintf("   cellchat_young.rds, cellchat_aged.rds, cellchat_merged.rds\n"))
cat(sprintf("   figures/ — overall comparison\n"))
cat(sprintf("   microglia/ — Microglia-focused\n"))
