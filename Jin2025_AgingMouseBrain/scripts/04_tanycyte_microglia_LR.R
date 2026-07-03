# =============================================================================
# 提取 Tanycyte ↔ Microglia 配体-受体互作对
# =============================================================================
suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
})

RESULT_DIR <- "d:/vscode/Jin2025_AgingMouseBrain/results/cellchat"

# Load CellChat objects
cat("Loading CellChat objects...\n")
cc_y <- readRDS(file.path(RESULT_DIR, "cellchat_young.rds"))
cc_a <- readRDS(file.path(RESULT_DIR, "cellchat_aged.rds"))

# Get cell type indices
ct_y <- levels(cc_y@idents)
ct_a <- levels(cc_a@idents)

cat("\nCell types:", paste(ct_y, collapse=", "), "\n")

tany_idx_y <- which(ct_y == "Tanycyte")
tany_idx_a <- which(ct_a == "Tanycyte")
mg_idx_y <- which(ct_y == "Microglia")
mg_idx_a <- which(ct_a == "Microglia")

# Extract LR pairs from the CellChat database
LR_db <- cc_y@DB$interaction
cat(sprintf("\nLR database: %d interactions\n", nrow(LR_db)))

# Get the inferred LR pairs from the CellChat objects
# cc@net$prob gives the communication probability matrix

# ---- Function to extract top LR pairs for a specific sender->receiver ----
extract_lr_pairs <- function(cc, sender_idx, receiver_idx, label) {
  cat(sprintf("\n=== %s ===\n", label))

  # Get the LR pair dataframe from the CellChat object
  # The inferred interactions are stored in a specific format

  # Try to access the inferred LR pair info
  if (length(cc@net) == 0 || !"prob" %in% names(cc@net)) {
    cat("  No probability matrix found\n")
    return(NULL)
  }

  # Get the communication probability for this specific sender-receiver pair
  prob <- cc@net$prob

  # The prob matrix is [sources, targets, LR_pairs]
  # We need to extract LR pairs for the specific sender-receiver

  # Get LR pair names
  lr_names <- dimnames(prob)[[3]]
  if (is.null(lr_names)) {
    cat("  No LR pair names in prob matrix\n")
    return(NULL)
  }

  # Extract probabilities for Tany->MG and MG->Tany
  tany_to_mg <- prob[tany_idx_y, mg_idx_y, ]
  mg_to_tany <- prob[mg_idx_y, tany_idx_y, ]

  # Create data frame of LR pairs sorted by probability
  t2m_df <- data.frame(
    LR_pair = lr_names,
    direction = "Tanycyte→Microglia",
    prob = as.numeric(tany_to_mg),
    stringsAsFactors = FALSE
  )
  t2m_df <- t2m_df[order(t2m_df$prob, decreasing = TRUE), ]

  m2t_df <- data.frame(
    LR_pair = lr_names,
    direction = "Microglia→Tanycyte",
    prob = as.numeric(mg_to_tany),
    stringsAsFactors = FALSE
  )
  m2t_df <- m2t_df[order(m2t_df$prob, decreasing = TRUE), ]

  # Filter to non-zero probabilities
  t2m_sig <- t2m_df[t2m_df$prob > 0, ]
  m2t_sig <- m2t_df[m2t_df$prob > 0, ]

  cat(sprintf("  Tany→MG: %d LR pairs with prob > 0\n", nrow(t2m_sig)))
  cat(sprintf("  MG→Tany: %d LR pairs with prob > 0\n", nrow(m2t_sig)))

  # Print top 20
  if (nrow(t2m_sig) > 0) {
    cat("\n  Top 20 Tanycyte → Microglia LR pairs:\n")
    top <- head(t2m_sig, 20)
    for (i in 1:nrow(top)) {
      cat(sprintf("    %2d. %-60s prob=%.2e\n", i, top$LR_pair[i], top$prob[i]))
    }
  }

  if (nrow(m2t_sig) > 0) {
    cat("\n  Top 20 Microglia → Tanycyte LR pairs:\n")
    top <- head(m2t_sig, 20)
    for (i in 1:nrow(top)) {
      cat(sprintf("    %2d. %-60s prob=%.2e\n", i, top$LR_pair[i], top$prob[i]))
    }
  }

  # Return combined dataframe
  result <- rbind(t2m_sig, m2t_sig)
  return(result)
}

# ---- Young ----
cat("\n\n========== YOUNG (2m) ==========")
young_lr <- extract_lr_pairs(cc_y, tany_idx_y, mg_idx_y, "Young(2m)")

# ---- Aged ----
cat("\n\n========== AGED (18m) ==========")
aged_lr <- extract_lr_pairs(cc_a, tany_idx_a, mg_idx_a, "Aged(18m)")

# ---- Compare: LR pairs that change with aging ----
cat("\n\n========== AGING CHANGES (Aged - Young) ==========\n")

if (!is.null(young_lr) && !is.null(aged_lr)) {
  # Merge young and aged
  young_lr$key <- paste(young_lr$LR_pair, young_lr$direction, sep="|")
  aged_lr$key <- paste(aged_lr$LR_pair, aged_lr$direction, sep="|")

  combined <- merge(young_lr[, c("key", "prob")], aged_lr[, c("key", "prob")],
                    by = "key", all = TRUE, suffixes = c("_young", "_aged"))
  combined[is.na(combined)] <- 0
  combined$change <- combined$prob_aged - combined$prob_young
  combined$log2FC <- log2((combined$prob_aged + 1e-30) / (combined$prob_young + 1e-30))

  # Split key back
  parts <- strsplit(combined$key, "\\|")
  combined$LR_pair <- sapply(parts, `[`, 1)
  combined$direction <- sapply(parts, `[`, 2)

  combined <- combined[order(abs(combined$change), decreasing = TRUE), ]

  # Save
  fwrite(combined[, c("LR_pair", "direction", "prob_young", "prob_aged", "change", "log2FC")],
         file.path(RESULT_DIR, "microglia", "tanycyte_microglia_LR_changes.csv"))
  cat("Saved: tanycyte_microglia_LR_changes.csv\n")

  # Show top changes
  cat("\nTop 20 increased LR pairs (Aged > Young):\n")
  top_up <- head(combined[combined$change > 0, ], 20)
  for (i in 1:nrow(top_up)) {
    cat(sprintf("  %2d. %-55s %-25s young=%.2e aged=%.2e log2FC=%+.2f\n",
                i, top_up$LR_pair[i], top_up$direction[i],
                top_up$prob_young[i], top_up$prob_aged[i], top_up$log2FC[i]))
  }

  cat("\nTop 20 decreased LR pairs (Aged < Young):\n")
  top_down <- head(combined[combined$change < 0, ], 20)
  for (i in 1:nrow(top_down)) {
    cat(sprintf("  %2d. %-55s %-25s young=%.2e aged=%.2e log2FC=%+.2f\n",
                i, top_down$LR_pair[i], top_down$direction[i],
                top_down$prob_young[i], top_down$prob_aged[i], top_down$log2FC[i]))
  }

  # Also save individual direction files
  for (dir in c("Tanycyte→Microglia", "Microglia→Tanycyte")) {
    dir_label <- if(dir == "Tanycyte→Microglia") "Tany_to_MG" else "MG_to_Tany"
    subset <- combined[combined$direction == dir, ]
    subset <- subset[order(abs(subset$change), decreasing = TRUE), ]
    fwrite(subset, file.path(RESULT_DIR, "microglia",
           paste0(dir_label, "_LR_changes.csv")))
  }
}

cat("\n\nDone!\n")
