## 00_install_packages.R
## Installs every package the rest of analysis_dge_host_pathogen/ depends on.
## Run once: source("00_install_packages.R"); manage_packages("automatic")

cran_packages <- c("ggplot2", "dplyr", "tidyr", "tibble", "pheatmap",
                    "RColorBrewer", "patchwork", "ggrepel", "data.table")
bioc_packages <- c("DESeq2")

manage_packages <- function(mode = c("interactive", "automatic")) {
  mode <- match.arg(mode)
  missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
  missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_cran) == 0 && length(missing_bioc) == 0) {
    message("All required packages are already installed.")
    return(invisible(NULL))
  }

  message("Missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  message("Missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))

  if (mode == "interactive") {
    ans <- readline("Install missing packages now? [y/N] ")
    if (!tolower(ans) %in% c("y", "yes")) {
      message("Skipped installation.")
      return(invisible(NULL))
    }
  }

  if (length(missing_cran) > 0) install.packages(missing_cran)
  if (length(missing_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
  }
}
