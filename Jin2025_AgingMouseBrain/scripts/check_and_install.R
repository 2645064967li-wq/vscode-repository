options(repos = c(CRAN = "https://cloud.r-project.org"))

# Check all CellChat deps
all_deps <- c("svglite", "expm", "sna", "ggpubr", "ggnetwork",
              "Biobase", "ComplexHeatmap", "NMF", "igraph", "Rcpp",
              "RcppArmadillo", "circlize", "ggalluvial", "reshape2",
              "dplyr", "ggplot2", "patchwork", "data.table")

for (pkg in all_deps) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("%-20s %s\n", pkg, if(ok) "OK" else "MISSING"))
}

# Install missing ones
for (pkg in c("svglite", "expm", "sna", "ggpubr", "ggnetwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    tryCatch({
      install.packages(pkg, type = "binary")
      cat("  Done\n")
    }, error = function(e) {
      cat("  FAILED:", e$message, "\n")
    })
  }
}

# Install CellChat from local source
if (!requireNamespace("CellChat", quietly = TRUE)) {
  cat("Installing CellChat from local source...\n")
  tryCatch({
    install.packages("d:/vscode/Jin2025_AgingMouseBrain/sqjin-CellChat-e4f6862",
                     repos = NULL, type = "source")
    cat("  Done\n")
  }, error = function(e) {
    cat("  FAILED:", e$message, "\n")
  })
}

cat("\n=== Final Check ===\n")
cat("CellChat installed:", requireNamespace("CellChat", quietly = TRUE), "\n")
