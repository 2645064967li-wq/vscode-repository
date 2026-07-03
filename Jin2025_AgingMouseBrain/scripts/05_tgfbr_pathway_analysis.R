# =============================================================================
# TGFBR通路深入分析: Tanycyte ↔ Microglia
# =============================================================================
# 重点分析 TGF-beta 配体-受体互作在衰老中的变化
# 包含: TGFB1/2/3 - TGFBR1/2/3, ACVR, BMPR 等
# =============================================================================
suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(data.table)
  library(reshape2)
})

RESULT_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/results/cellchat"
OUT_DIR <- file.path(RESULT_DIR, "tgfbr_pathway")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Loading ----
cat("Loading CellChat objects...\n")
cc_y <- readRDS(file.path(RESULT_DIR, "cellchat_young.rds"))
cc_a <- readRDS(file.path(RESULT_DIR, "cellchat_aged.rds"))

ct_y <- levels(cc_y@idents)
ct_a <- levels(cc_a@idents)
cat("Cell types:", paste(ct_y, collapse=", "), "\n")

# ---- 2. Extract ALL TGF-beta superfamily interactions ----
cat("\n=== Extracting TGF-beta superfamily LR pairs ===\n")

# The CellChat database has pathway annotations
# TGFb signaling pathway includes TGFB1/2/3, BMPs, Activins, etc.
db <- cc_y@DB$interaction

# Find all interactions in TGFb pathway
tgf_pathways <- unique(db$pathway_name[grep("TGF|TGFb|BMP|ACTIVIN|NODAL",
                                             db$pathway_name, ignore.case=TRUE)])
cat("TGFb-related pathways in DB:", paste(tgf_pathways, collapse=", "), "\n")

# Get all LR pairs annotated as TGFb
tgf_lr <- db[db$pathway_name %in% tgf_pathways, ]
cat(sprintf("Total TGFb-related LR pairs in DB: %d\n", nrow(tgf_lr)))

# Also find any LR pairs whose ligand or receptor contains TGF/ACVR/BMPR
tgf_related <- db[grepl("TGF|TGFB|ACVR|BMPR|AMHR|ENG|TDGF|GDF|BMP|INHB|NODAL|L EF|TGFBR",
                         db$interaction_name, ignore.case=TRUE), ]
cat(sprintf("Extended TGF-related LR pairs: %d\n", nrow(tgf_related)))

# ---- 3. Extract Tany-MG TGFBR interactions from CellChat objects ----
cat("\n=== Tanycyte <-> Microglia TGFBR interactions ===\n")

tany_idx <- which(ct_y == "Tanycyte")
mg_idx <- which(ct_y == "Microglia")

extract_tgfb_lr <- function(cc, label) {
  prob <- cc@net$prob
  lr_names <- dimnames(prob)[[3]]

  # Find TGF-related LR pairs
  tgf_pattern <- "TGF|TGFB|ACVR|BMPR|ENG|TDGF|GDF|BMP[0-9]|INHB|NODAL|LEFTY|AMH|TGFR"
  tgf_idx <- grep(tgf_pattern, lr_names, ignore.case=TRUE)

  if (length(tgf_idx) == 0) {
    cat(sprintf("  [%s] No TGF-related LR pairs found\n", label))
    return(NULL)
  }

  cat(sprintf("  [%s] %d TGF-related LR pairs found in CellChat object\n", label, length(tgf_idx)))

  # Extract probabilities
  result <- data.frame(
    LR_pair = lr_names[tgf_idx],
    Tany_to_MG = as.numeric(prob[tany_idx, mg_idx, tgf_idx]),
    MG_to_Tany = as.numeric(prob[mg_idx, tany_idx, tgf_idx]),
    stringsAsFactors = FALSE
  )

  # Filter to those with at least one non-zero probability
  result <- result[result$Tany_to_MG > 0 | result$MG_to_Tany > 0, ]
  result <- result[order(pmax(result$Tany_to_MG, result$MG_to_Tany), decreasing=TRUE), ]

  cat(sprintf("  Non-zero Tany<->MG interactions: %d\n", nrow(result)))

  if (nrow(result) > 0) {
    cat("\n  Top TGFBR LR pairs:\n")
    for (i in 1:min(20, nrow(result))) {
      cat(sprintf("    %2d. %-50s Tany->MG=%.2e  MG->Tany=%.2e\n",
                  i, result$LR_pair[i], result$Tany_to_MG[i], result$MG_to_Tany[i]))
    }
  }

  return(result)
}

young_tgf <- extract_tgfb_lr(cc_y, "Young(2m)")
aged_tgf <- extract_tgfb_lr(cc_a, "Aged(18m)")

# ---- 4. Differential analysis ----
cat("\n=== Differential TGFBR signaling (Aged - Young) ===\n")

if (!is.null(young_tgf) && !is.null(aged_tgf)) {
  # Merge
  merged <- merge(young_tgf, aged_tgf, by="LR_pair", all=TRUE,
                  suffixes=c("_young", "_aged"))
  merged[is.na(merged)] <- 0

  # Calculate changes
  merged$Tany_to_MG_change <- merged$Tany_to_MG_aged - merged$Tany_to_MG_young
  merged$MG_to_Tany_change <- merged$MG_to_Tany_aged - merged$MG_to_Tany_young
  merged$max_prob <- pmax(merged$Tany_to_MG_young, merged$Tany_to_MG_aged,
                          merged$MG_to_Tany_young, merged$MG_to_Tany_aged)

  # Log2FC
  eps <- 1e-30
  merged$Tany_to_MG_log2FC <- log2((merged$Tany_to_MG_aged + eps) /
                                     (merged$Tany_to_MG_young + eps))
  merged$MG_to_Tany_log2FC <- log2((merged$MG_to_Tany_aged + eps) /
                                     (merged$MG_to_Tany_young + eps))

  # Overall change score
  merged$total_change <- abs(merged$Tany_to_MG_change) + abs(merged$MG_to_Tany_change)

  merged <- merged[order(merged$total_change, decreasing=TRUE), ]

  # Save
  fwrite(merged, file.path(OUT_DIR, "tgfbr_Tany_MG_differential.csv"))

  # ---- Print results ----
  cat("\n--- TGFBR: Increased with aging ---\n")
  up <- merged[merged$Tany_to_MG_change > 0 | merged$MG_to_Tany_change > 0, ]
  if (nrow(up) > 0) {
    for (i in 1:min(15, nrow(up))) {
      cat(sprintf("  %2d. %-45s Tany->MG: %.2e->%.2e (FC=%+.2f)  MG->Tany: %.2e->%.2e (FC=%+.2f)\n",
                  i, up$LR_pair[i],
                  up$Tany_to_MG_young[i], up$Tany_to_MG_aged[i], up$Tany_to_MG_log2FC[i],
                  up$MG_to_Tany_young[i], up$MG_to_Tany_aged[i], up$MG_to_Tany_log2FC[i]))
    }
  }

  cat("\n--- TGFBR: Decreased with aging ---\n")
  down <- merged[merged$Tany_to_MG_change < 0 | merged$MG_to_Tany_change < 0, ]
  if (nrow(down) > 0) {
    for (i in 1:min(15, nrow(down))) {
      cat(sprintf("  %2d. %-45s Tany->MG: %.2e->%.2e (FC=%+.2f)  MG->Tany: %.2e->%.2e (FC=%+.2f)\n",
                  i, down$LR_pair[i],
                  down$Tany_to_MG_young[i], down$Tany_to_MG_aged[i], down$Tany_to_MG_log2FC[i],
                  down$MG_to_Tany_young[i], down$MG_to_Tany_aged[i], down$MG_to_Tany_log2FC[i]))
    }
  }
}

# ---- 5. Extract specific TGFBR ligand-receptor sub-components ----
cat("\n=== Key TGFBR Ligand-Receptor Breakdown ===\n")

# The CellChat LR names encode: Ligand_Receptor or Ligand_Receptor1_Receptor2
# Key TGF-beta components:
# Ligands: TGFB1, TGFB2, TGFB3, GDF11, INHBA, INHBB, BMPs
# Receptors: TGFBR1, TGFBR2, TGFBR3, ACVR1, ACVR2A/B, BMPR1A/B, BMPR2

key_ligands <- c("TGFB1", "TGFB2", "TGFB3", "GDF11", "INHBA", "INHBB",
                 "BMP2", "BMP4", "BMP6", "BMP7", "NODAL", "LEFTY1", "LEFTY2")
key_receptors <- c("TGFBR1", "TGFBR2", "TGFBR3", "ACVR1", "ACVR1B", "ACVR1C",
                   "ACVR2A", "ACVR2B", "BMPR1A", "BMPR1B", "BMPR2", "ENG", "TDGF1")

cat("\nKey TGF-beta Ligands:", paste(key_ligands, collapse=", "))
cat("\nKey TGF-beta Receptors:", paste(key_receptors, collapse=", "))

# ---- 6. NetVisual: TGFb pathway specific ----
cat("\n\n=== Generating TGFBR pathway visualizations ===\n")

# Try to use CellChat's signaling pathway visualization
# First, identify the TGFb pathway index
# Note: computeCommunProbPathway was buggy, so we work at LR-level

# Manual approach: aggregate TGFBR-specific interactions
if (!is.null(merged)) {

  # ---- 6a. Create summary bar plot ----
  cat("Creating TGFBR summary plots...\n")

  # Prepare data for visualization
  plot_data <- merged
  plot_data$LR_short <- gsub("_.*", "", plot_data$LR_pair)  # Get ligand name

  # Top 15 TGFBR interactions by max probability
  top_15 <- head(plot_data[order(plot_data$max_prob, decreasing=TRUE), ], 15)

  # Bar plot comparing young vs aged
  png(file.path(OUT_DIR, "01_tgfbr_TanyMG_barplot.png"),
      width=1600, height=900, res=150)

  par(mar=c(5, 18, 4, 2))

  # Tany -> MG
  layout(matrix(1:2, 1, 2))

  # Panel A: Tanycyte -> Microglia
  bp_data <- top_15[order(top_15$Tany_to_MG_aged), ]
  y_pos <- seq_len(nrow(bp_data))

  plot(bp_data$Tany_to_MG_young, y_pos, type="p", pch=19, col="#4DBBD5", cex=1.5,
       xlim=range(0, max(bp_data$Tany_to_MG_young, bp_data$Tany_to_MG_aged)*1.2),
       ylim=c(0.5, nrow(bp_data)+0.5), yaxt="n",
       xlab="Communication Probability", ylab="",
       main="Tanycyte -> Microglia (TGFBR)")
  points(bp_data$Tany_to_MG_aged, y_pos, pch=19, col="#E64B35", cex=1.5)
  axis(2, at=y_pos, labels=bp_data$LR_pair, las=2, cex.axis=0.7)
  legend("bottomright", c("Young(2m)", "Aged(18m)"),
         col=c("#4DBBD5", "#E64B35"), pch=19, cex=0.9)

  # Panel B: Microglia -> Tanycyte
  bp_data2 <- top_15[order(top_15$MG_to_Tany_aged), ]

  plot(bp_data2$MG_to_Tany_young, y_pos, type="p", pch=19, col="#4DBBD5", cex=1.5,
       xlim=range(0, max(bp_data2$MG_to_Tany_young, bp_data2$MG_to_Tany_aged)*1.2),
       ylim=c(0.5, nrow(bp_data2)+0.5), yaxt="n",
       xlab="Communication Probability", ylab="",
       main="Microglia -> Tanycyte (TGFBR)")
  points(bp_data2$MG_to_Tany_aged, y_pos, pch=19, col="#E64B35", cex=1.5)
  axis(2, at=y_pos, labels=bp_data2$LR_pair, las=2, cex.axis=0.7)
  legend("bottomright", c("Young(2m)", "Aged(18m)"),
         col=c("#4DBBD5", "#E64B35"), pch=19, cex=0.9)

  dev.off()
  cat("  [OK] 01_tgfbr_TanyMG_barplot.png\n")

  # ---- 6b. Volcano-style plot (log2FC vs max prob) ----
  png(file.path(OUT_DIR, "02_tgfbr_volcano.png"),
      width=1400, height=1000, res=150)

  par(mfrow=c(1, 2), mar=c(5,5,4,2))

  # Tany->MG volcano
  plot(merged$Tany_to_MG_log2FC, -log10(merged$max_prob + 1e-30),
       pch=20, col=rgb(0.5, 0.5, 0.5, 0.5), cex=1,
       xlab="log2FC (Aged/Young)", ylab="-log10(max prob)",
       main="Tanycyte -> Microglia: TGFBR Changes")
  abline(v=0, lty=2, col="grey")
  abline(h=-log10(0.05), lty=3, col="red")

  # Highlight significant changes
  sig_up <- merged[merged$Tany_to_MG_log2FC > 0.5, ]
  sig_down <- merged[merged$Tany_to_MG_log2FC < -0.5, ]

  if (nrow(sig_up) > 0) {
    points(sig_up$Tany_to_MG_log2FC, -log10(sig_up$max_prob + 1e-30),
           pch=20, col="#E64B35", cex=1.3)
    text(sig_up$Tany_to_MG_log2FC, -log10(sig_up$max_prob + 1e-30),
         labels=sig_up$LR_pair, pos=4, cex=0.5, col="#E64B35")
  }
  if (nrow(sig_down) > 0) {
    points(sig_down$Tany_to_MG_log2FC, -log10(sig_down$max_prob + 1e-30),
           pch=20, col="#4DBBD5", cex=1.3)
    text(sig_down$Tany_to_MG_log2FC, -log10(sig_down$max_prob + 1e-30),
         labels=sig_down$LR_pair, pos=2, cex=0.5, col="#4DBBD5")
  }
  legend("topright", c("Up in Aged", "Down in Aged", "No change"),
         col=c("#E64B35", "#4DBBD5", "grey"), pch=20, cex=0.8)

  # MG->Tany volcano
  plot(merged$MG_to_Tany_log2FC, -log10(merged$max_prob + 1e-30),
       pch=20, col=rgb(0.5, 0.5, 0.5, 0.5), cex=1,
       xlab="log2FC (Aged/Young)", ylab="-log10(max prob)",
       main="Microglia -> Tanycyte: TGFBR Changes")
  abline(v=0, lty=2, col="grey")
  abline(h=-log10(0.05), lty=3, col="red")

  sig_up2 <- merged[merged$MG_to_Tany_log2FC > 0.5, ]
  sig_down2 <- merged[merged$MG_to_Tany_log2FC < -0.5, ]

  if (nrow(sig_up2) > 0) {
    points(sig_up2$MG_to_Tany_log2FC, -log10(sig_up2$max_prob + 1e-30),
           pch=20, col="#E64B35", cex=1.3)
    text(sig_up2$MG_to_Tany_log2FC, -log10(sig_up2$max_prob + 1e-30),
         labels=sig_up2$LR_pair, pos=4, cex=0.5, col="#E64B35")
  }
  if (nrow(sig_down2) > 0) {
    points(sig_down2$MG_to_Tany_log2FC, -log10(sig_down2$max_prob + 1e-30),
           pch=20, col="#4DBBD5", cex=1.3)
    text(sig_down2$MG_to_Tany_log2FC, -log10(sig_down2$max_prob + 1e-30),
         labels=sig_down2$LR_pair, pos=2, cex=0.5, col="#4DBBD5")
  }
  legend("topright", c("Up in Aged", "Down in Aged", "No change"),
         col=c("#E64B35", "#4DBBD5", "grey"), pch=20, cex=0.8)

  dev.off()
  cat("  [OK] 02_tgfbr_volcano.png\n")

  # ---- 6c. TGFBR-specific LR pair heatmap data ----
  # Create a matrix of all TGFBR interactions across cell types
  cat("\nCreating TGFBR cross-cell-type summary...\n")
}

# ---- 7. TGF-beta ligands/receptors: global pattern across all cell types ----
cat("\n=== TGFBR signaling: Global cross-cell-type pattern ===\n")

# Extract TGFBR interactions across ALL cell type pairs, not just Tany-MG
extract_tgfb_all_pairs <- function(cc, label) {
  prob <- cc@net$prob
  lr_names <- dimnames(prob)[[3]]
  tgf_pattern <- "TGF|TGFB|ACVR|BMPR|ENG|TDGF|GDF|BMP[0-9]|INHB|NODAL|LEFTY|AMH|TGFR"
  tgf_idx <- grep(tgf_pattern, lr_names, ignore.case=TRUE)

  if (length(tgf_idx) == 0) return(NULL)

  # Subset to TGFBR interactions
  tgf_prob <- prob[,, tgf_idx, drop=FALSE]

  # Sum across all TGFBR LRs for each cell type pair
  n_ct <- dim(prob)[1]
  ct_names <- dimnames(prob)[[1]]

  # Total TGFBR signaling strength per cell pair
  tgf_sum <- apply(tgf_prob, c(1, 2), sum)
  rownames(tgf_sum) <- ct_names
  colnames(tgf_sum) <- ct_names

  cat(sprintf("\n  [%s] TGFBR total signaling:\n", label))
  cat(sprintf("    Max sender: %s\n", ct_names[which.max(rowSums(tgf_sum))]))
  cat(sprintf("    Max receiver: %s\n", ct_names[which.max(colSums(tgf_sum))]))

  # Tanycyte-specific
  tany_send <- sum(tgf_sum["Tanycyte", ])
  tany_recv <- sum(tgf_sum[, "Tanycyte"])
  cat(sprintf("    Tanycyte outgoing: %.2e\n", tany_send))
  cat(sprintf("    Tanycyte incoming: %.2e\n", tany_recv))

  return(tgf_sum)
}

tgf_y_all <- extract_tgfb_all_pairs(cc_y, "Young(2m)")
tgf_a_all <- extract_tgfb_all_pairs(cc_a, "Aged(18m)")

if (!is.null(tgf_y_all) && !is.null(tgf_a_all)) {
  # Differential TGFBR matrix
  tgf_diff <- tgf_a_all - tgf_y_all

  # Save
  write.csv(tgf_y_all, file.path(OUT_DIR, "tgfbr_matrix_young.csv"))
  write.csv(tgf_a_all, file.path(OUT_DIR, "tgfbr_matrix_aged.csv"))
  write.csv(tgf_diff, file.path(OUT_DIR, "tgfbr_matrix_diff.csv"))

  # Print Tany-related changes
  tany_ct <- "Tanycyte"
  mg_ct <- "Microglia"

  cat("\n--- TGFBR: Tanycyte outgoing changes (all targets) ---\n")
  tany_out <- tgf_diff[tany_ct, ]
  tany_out <- sort(tany_out, decreasing=TRUE)
  for (i in seq_along(tany_out)) {
    if (abs(tany_out[i]) > 0) {
      direction <- ifelse(tany_out[i] > 0, "+", "")
      cat(sprintf("  Tany->%-25s %s%.2e\n", names(tany_out)[i], direction, tany_out[i]))
    }
  }

  cat("\n--- TGFBR: Microglia incoming changes (all sources) ---\n")
  mg_in <- tgf_diff[, mg_ct]
  mg_in <- sort(mg_in, decreasing=TRUE)
  for (i in seq_along(mg_in)) {
    if (abs(mg_in[i]) > 0) {
      direction <- ifelse(mg_in[i] > 0, "+", "")
      cat(sprintf("  %-25s ->MG  %s%.2e\n", names(mg_in)[i], direction, mg_in[i]))
    }
  }
}

# ---- 8. Summary Report ----
cat("\n\n========================================\n")
cat(" TGFBR Pathway Analysis Complete\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", OUT_DIR))
cat("Files:\n")
cat("  - tgfbr_Tany_MG_differential.csv\n")
cat("  - tgfbr_matrix_young.csv / aged.csv / diff.csv\n")
cat("  - 01_tgfbr_TanyMG_barplot.png\n")
cat("  - 02_tgfbr_volcano.png\n")
