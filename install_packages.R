cran_packages <- c(
  "shiny",
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "circlize",
  "openxlsx",
  "readxl",
  "DT",
  "tibble",
  "scales"
)

install.packages(cran_packages, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
