# Install CellChat dependencies
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install CRAN packages
cran_pkgs <- c("dplyr", "ggplot2", "patchwork", "data.table", "Matrix",
               "circlize", "RColorBrewer", "ggalluvial", "igraph", "NMF",
               "reshape2", "tidyr", "plyr", "Rcpp", "RcppArmadillo",
               "foreach", "doParallel", "Rtsne", "scales", "magrittr")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, type = "binary")
  } else {
    cat(pkg, "already installed\n")
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_pkgs <- c("Biobase", "BiocGenerics", "ComplexHeatmap")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(pkg, "already installed\n")
  }
}

# Install CellChat from GitHub
if (!requireNamespace("CellChat", quietly = TRUE)) {
  cat("Installing CellChat from GitHub...\n")
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("sqjin/CellChat", upgrade = "never")
} else {
  cat("CellChat already installed\n")
}

cat("\nAll CellChat dependencies installed!\n")
