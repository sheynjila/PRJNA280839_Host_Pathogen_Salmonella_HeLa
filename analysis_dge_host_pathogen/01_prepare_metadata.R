## 01_prepare_metadata.R
## Loads the independently-verified GEO/SRA metadata for this project's
## real 8-sample subset (WT Salmonella infection vs mock, 4h and 24h,
## 2 replicates each -- all single-end, verified against SRA run metadata).

metadata <- read.csv("../data/metadata_8samples.csv", stringsAsFactors = FALSE)

cat("Sample metadata (independently verified against SRA run info):\n")
print(metadata[, c("gsm", "condition", "timepoint_h", "replicate", "published_column")])

metadata$condition <- factor(metadata$condition, levels = c("mock", "infected"))
metadata$timepoint_h <- factor(metadata$timepoint_h, levels = c(4, 24))
rownames(metadata) <- metadata$published_column

cat("\nDesign: 2 conditions (mock/infected) x 2 timepoints (4h/24h) x 2 replicates = 8 samples\n")
