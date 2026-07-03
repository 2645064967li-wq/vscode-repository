# =============================================================================
# CellChat Analysis v2: 跳过computeCommunProbPathway，直接比较LR-level网络
# =============================================================================
set.seed(42)
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("\n============================================\n")
cat(" CellChat v2: LR-level comparison\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
})

DATA_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/data/cellchat"
RESULT_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/results/cellchat"
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULT_DIR, "microglia"), recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load data ----
cat("\n[1] Loading data...\n")
young_meta <- fread(file.path(DATA_DIR, "young_metadata.csv"))
aged_meta  <- fread(file.path(DATA_DIR, "aged_metadata.csv"))

# Filter Unclassified & small types
MIN_CELLS <- 20
young_meta <- young_meta[young_meta$cell_type != "Unclassified", ]
aged_meta  <- aged_meta[aged_meta$cell_type != "Unclassified", ]
ct_y <- young_meta[, .N, by = cell_type][N >= MIN_CELLS, cell_type]
ct_a <- aged_meta[, .N, by = cell_type][N >= MIN_CELLS, cell_type]
keep_ct <- intersect(ct_y, ct_a)
young_meta <- young_meta[young_meta$cell_type %in% keep_ct, ]
aged_meta  <- aged_meta[aged_meta$cell_type %in% keep_ct, ]
cat(sprintf("Young: %d cells, %d types\n", nrow(young_meta), length(keep_ct)))
cat(sprintf("Aged:  %d cells, %d types\n", nrow(aged_meta), length(keep_ct)))

# ---- 2. Build CellChat objects ----
cat("\n[2] Creating CellChat objects...\n")

build_cc <- function(label, meta) {
  cat(sprintf("\n  --- %s ---\n", label))
  counts_file <- file.path(DATA_DIR, paste0(label, "_raw_counts.csv.gz"))
  cat("    Loading counts...\n")
  counts_dt <- fread(counts_file)
  gene_names <- counts_dt[[1]]
  counts_mat <- as.matrix(counts_dt[, -1, drop = FALSE])
  rownames(counts_mat) <- gene_names

  valid_cells <- intersect(colnames(counts_mat), meta$cell_label)
  counts_mat <- counts_mat[, valid_cells, drop = FALSE]
  meta_use <- meta[match(valid_cells, meta$cell_label), ]

  cat(sprintf("    %d genes x %d cells\n", nrow(counts_mat), ncol(counts_mat)))

  cc <- createCellChat(object = counts_mat, meta = meta_use, group.by = "cell_type")
  cc@DB <- CellChatDB.mouse
  cc <- subsetData(cc)
  cat(sprintf("    %d signaling genes\n", nrow(cc@data.signaling)))
  rm(counts_dt, counts_mat)
  gc()
  return(cc)
}

cc_y <- build_cc("young", young_meta)
cc_a <- build_cc("aged", aged_meta)

# ---- 3. Run core CellChat pipeline (no pathway) ----
cat("\n[3] Running CellChat pipeline...\n")

run_chat <- function(cc, label) {
  cat(sprintf("\n  === %s ===\n", label))
  cat("    identifyOverExpressedGenes...\n")
  cc <- identifyOverExpressedGenes(cc)
  cat("    identifyOverExpressedInteractions...\n")
  cc <- identifyOverExpressedInteractions(cc)
  cat("    computeCommunProb...\n")
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                           raw.use = TRUE, population.size = TRUE)
  cat("    filterCommunication...\n")
  cc <- filterCommunication(cc, min.cells = 10)
  cat("    aggregateNet...\n")
  cc <- aggregateNet(cc)
  return(cc)
}

cc_y <- run_chat(cc_y, "Young(2m)")
cc_a <- run_chat(cc_a, "Aged(18m)")

# Save
saveRDS(cc_y, file.path(RESULT_DIR, "cellchat_young.rds"))
saveRDS(cc_a, file.path(RESULT_DIR, "cellchat_aged.rds"))
cat("\n  CellChat objects saved\n")

# ---- 4. Compare networks at LR interaction level ----
cat("\n[4] Comparing networks...\n")

# Get common cell types
common_ct <- intersect(levels(cc_y@idents), levels(cc_a@idents))
cat(sprintf("  Common cell types: %d\n", length(common_ct)))

# Net interaction counts
net_y <- cc_y@net$count
net_a <- cc_a@net$count

cat("\n  === Total interactions ===\n")
cat(sprintf("  Young: %.0f interactions\n", sum(net_y)))
cat(sprintf("  Aged:  %.0f interactions\n", sum(net_a)))
cat(sprintf("  Ratio (Aged/Young): %.2f\n", sum(net_a) / sum(net_y)))

# Per cell type interaction export
write_net_stats <- function(cc, net, label) {
  ct <- levels(cc@idents)
  df <- data.frame(
    cell_type = ct,
    outgoing_count = colSums(as.matrix(net)),
    incoming_count = rowSums(as.matrix(net)),
    outgoing_weight = colSums(as.matrix(cc@net$weight)),
    incoming_weight = rowSums(as.matrix(cc@net$weight))
  )
  df <- df[order(df$outgoing_count, decreasing = TRUE), ]
  fwrite(df, file.path(RESULT_DIR, paste0(label, "_net_stats.csv")))
  cat(sprintf("    %s_net_stats.csv saved\n", label))
  return(df)
}

cat("\n  Young per-type stats:\n")
young_stats <- write_net_stats(cc_y, net_y, "young")
print(young_stats, row.names = FALSE)

cat("\n  Aged per-type stats:\n")
aged_stats <- write_net_stats(cc_a, net_a, "aged")
print(aged_stats, row.names = FALSE)

# Differential: aged - young per cell type
cat("\n  === Differential outgoing (Aged - Young) ===\n")
diff_stats <- merge(young_stats, aged_stats, by = "cell_type", suffixes = c("_y", "_a"))
diff_stats$outgoing_diff <- diff_stats$outgoing_count_a - diff_stats$outgoing_count_y
diff_stats$incoming_diff <- diff_stats$incoming_count_a - diff_stats$incoming_count_y
diff_stats <- diff_stats[order(abs(diff_stats$outgoing_diff), decreasing = TRUE), ]
print(diff_stats[, c("cell_type", "outgoing_count_y", "outgoing_count_a",
                      "outgoing_diff", "incoming_diff")], row.names = FALSE)
fwrite(diff_stats, file.path(RESULT_DIR, "differential_net_stats.csv"))

# ---- 5. Figures ----
cat("\n[5] Generating figures...\n")
fig_dir <- file.path(RESULT_DIR, "figures")

# 5a. Circle plot: Young vs Aged overall network
cat("  5a. Circle plots...\n")
pdf(file.path(fig_dir, "01_circle_young.pdf"), width = 10, height = 10)
netVisual_circle(net_y, weight.scale = TRUE,
                 title.name = "Young(2m) Interaction Network",
                 vertex.label.cex = 0.7)
dev.off()

pdf(file.path(fig_dir, "01_circle_aged.pdf"), width = 10, height = 10)
netVisual_circle(net_a, weight.scale = TRUE,
                 title.name = "Aged(18m) Interaction Network",
                 vertex.label.cex = 0.7)
dev.off()
cat("    OK\n")

# 5b. Heatmap comparison
cat("  5b. Heatmap...\n")
pdf(file.path(fig_dir, "02_heatmap_young.pdf"), width = 10, height = 8)
netVisual_heatmap(cc_y, measure = "count", font.size = 7,
                  title.name = "Young(2m): Interaction Counts")
dev.off()

pdf(file.path(fig_dir, "02_heatmap_aged.pdf"), width = 10, height = 8)
netVisual_heatmap(cc_a, measure = "count", font.size = 7,
                  title.name = "Aged(18m): Interaction Counts")
dev.off()

# Difference heatmap - skip due to CellChat API limitations
cat("    Diff heatmap skipped (CellChat API limitation)\n")

# 5c. Signaling role analysis (requires pathway data - may fail)
cat("  5c. Signaling role scatter...\n")
tryCatch({
  cc_y <- netAnalysis_computeCentrality(cc_y, slot.name = "netP")
  cc_a <- netAnalysis_computeCentrality(cc_a, slot.name = "netP")
  pdf(file.path(fig_dir, "03_signaling_role.pdf"), width = 14, height = 7)
  par(mfrow = c(1, 2))
  netAnalysis_signalingRole_scatter(cc_y, title = "Young(2m)", font.size = 8)
  netAnalysis_signalingRole_scatter(cc_a, title = "Aged(18m)", font.size = 8)
  dev.off()
  cat("    OK\n")
}, error = function(e) {
  cat("    Skipped (pathway not available):", conditionMessage(e), "\n")
})

# ---- 6. Microglia focus ----
cat("\n[6] Microglia-focused analysis...\n")
mg_dir <- file.path(RESULT_DIR, "microglia")

if ("Microglia" %in% common_ct) {
  other_ct <- setdiff(common_ct, "Microglia")
  mg_idx_y <- which(levels(cc_y@idents) == "Microglia")
  mg_idx_a <- which(levels(cc_a@idents) == "Microglia")
  other_y <- which(levels(cc_y@idents) %in% other_ct)
  other_a <- which(levels(cc_a@idents) %in% other_ct)

  # 6a. MG outgoing circle
  cat("  6a. MG outgoing...\n")
  tryCatch({
    pdf(file.path(mg_dir, "MG01_outgoing.pdf"), width = 14, height = 7)
    par(mfrow = c(1, 2))
    netVisual_circle(net_y, sources.use = mg_idx_y, targets.use = other_y,
      title.name = "Young: MG -> Others", vertex.label.cex = 0.7, remove.isolate = TRUE)
    netVisual_circle(net_a, sources.use = mg_idx_a, targets.use = other_a,
      title.name = "Aged: MG -> Others", vertex.label.cex = 0.7, remove.isolate = TRUE)
    dev.off()
    cat("    OK\n")
  }, error = function(e) {
    cat("    Skipped:", conditionMessage(e), "\n")
  })

  # 6b. MG incoming circle
  cat("  6b. MG incoming...\n")
  tryCatch({
    pdf(file.path(mg_dir, "MG02_incoming.pdf"), width = 14, height = 7)
    par(mfrow = c(1, 2))
    netVisual_circle(net_y, sources.use = other_y, targets.use = mg_idx_y,
      title.name = "Young: Others -> MG", vertex.label.cex = 0.7, remove.isolate = TRUE)
    netVisual_circle(net_a, sources.use = other_a, targets.use = mg_idx_a,
      title.name = "Aged: Others -> MG", vertex.label.cex = 0.7, remove.isolate = TRUE)
    dev.off()
    cat("    OK\n")
  }, error = function(e) {
    cat("    Skipped:", conditionMessage(e), "\n")
  })

  # 6c. MG outgoing targets ranked
  cat("  6c. MG target ranking...\n")
  mg_out_y <- as.matrix(net_y)[mg_idx_y, other_y]
  mg_out_a <- as.matrix(net_a)[mg_idx_a, other_a]
  names(mg_out_y) <- other_ct
  names(mg_out_a) <- other_ct

  mg_targets <- data.frame(
    target = other_ct,
    young_interactions = mg_out_y[other_ct],
    aged_interactions = mg_out_a[other_ct],
    change = mg_out_a[other_ct] - mg_out_y[other_ct]
  )
  mg_targets <- mg_targets[order(abs(mg_targets$change), decreasing = TRUE), ]
  fwrite(mg_targets, file.path(mg_dir, "MG03_microglia_targets.csv"))

  cat("\n  MG outgoing interaction changes (Aged - Young):\n")
  print(mg_targets, row.names = FALSE)

  # 6d. MG incoming from sources ranked
  cat("\n  6d. MG source ranking...\n")
  mg_in_y <- as.matrix(net_y)[other_y, mg_idx_y]
  mg_in_a <- as.matrix(net_a)[other_a, mg_idx_a]
  names(mg_in_y) <- other_ct
  names(mg_in_a) <- other_ct

  mg_sources <- data.frame(
    source = other_ct,
    young_interactions = mg_in_y[other_ct],
    aged_interactions = mg_in_a[other_ct],
    change = mg_in_a[other_ct] - mg_in_y[other_ct]
  )
  mg_sources <- mg_sources[order(abs(mg_sources$change), decreasing = TRUE), ]
  fwrite(mg_sources, file.path(mg_dir, "MG04_microglia_sources.csv"))

  cat("\n  MG incoming interaction changes (Aged - Young):\n")
  print(mg_sources, row.names = FALSE)

  # 6e. Chord diagram
  cat("  6e. Chord diagrams...\n")
  tryCatch({
    pdf(file.path(mg_dir, "MG05_chord.pdf"), width = 14, height = 7)
    par(mfrow = c(1, 2))
    netVisual_chord_gene(cc_y, sources.use = mg_idx_y, targets.use = other_y,
      title.name = "Young: MG Outgoing Signaling", legend.pos.x = 8)
    netVisual_chord_gene(cc_a, sources.use = mg_idx_a, targets.use = other_a,
      title.name = "Aged: MG Outgoing Signaling", legend.pos.x = 8)
    dev.off()
    cat("    OK\n")
  }, error = function(e) {
    cat("    Skipped:", conditionMessage(e), "\n")
  })

} else {
  cat("  Microglia NOT found in common cell types!\n")
}

# ---- 7. Summary ----
cat("\n============================================\n")
cat(" CellChat v2 analysis complete!\n")
cat("============================================\n")
cat(sprintf("\nYoung: %.0f total interactions\n", sum(net_y)))
cat(sprintf("Aged:  %.0f total interactions\n", sum(net_a)))
cat(sprintf("Difference: %.0f\n", sum(net_a) - sum(net_y)))

cat(sprintf("\nResults: %s/\n", RESULT_DIR))
cat("  cellchat_young.rds, cellchat_aged.rds\n")
cat("  young_net_stats.csv, aged_net_stats.csv\n")
cat("  differential_net_stats.csv\n")
cat("  figures/ — network visualization\n")
cat("  microglia/ — MG-focused analysis\n")
