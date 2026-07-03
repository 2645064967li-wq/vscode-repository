library(data.table)
cat("Testing data loading...\n")

# Test metadata
ymeta <- fread("d:/vscode/Jin2025_AgingMouseBrain/data/cellchat/young_metadata.csv")
cat(sprintf("Young metadata: %d cells, %d cell types\n",
            nrow(ymeta), length(unique(ymeta$cell_type))))

ameta <- fread("d:/vscode/Jin2025_AgingMouseBrain/data/cellchat/aged_metadata.csv")
cat(sprintf("Aged metadata: %d cells, %d cell types\n",
            nrow(ameta), length(unique(ameta$cell_type))))

# Test counts (read first few rows)
cat("Loading young counts header...\n")
ycounts <- fread("d:/vscode/Jin2025_AgingMouseBrain/data/cellchat/young_raw_counts.csv.gz", nrows = 5)
cat(sprintf("Young counts: %d cols (first col = gene names)\n", ncol(ycounts)))
cat("First 3 gene names:", head(ycounts[[1]], 3), "\n")

# Check if cell barcodes match between counts and metadata
y_cols <- colnames(ycounts)
y_cells_in_counts <- y_cols[-1]  # skip gene name col
y_cells_in_meta <- ymeta$cell_label
y_match <- sum(y_cells_in_counts %in% y_cells_in_meta)
cat(sprintf("Young: %d / %d cells in counts match metadata\n",
            y_match, length(y_cells_in_counts)))

cat("\nTest passed!\n")
