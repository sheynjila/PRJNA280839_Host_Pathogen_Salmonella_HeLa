## 02_load_and_validate_counts.R
## Loads the real published combined host+pathogen count matrix
## (GSE67758_raw_gene_quantifications.csv.gz) and splits it into a host
## count matrix and a pathogen count matrix, subset to this project's 8
## real samples.
##
## REAL FILE-FORMAT QUIRK, DISCLOSED: despite its ".csv" extension, this
## file is actually TAB-delimited, GTF/GFF-attribute-style (one row per
## genomic feature, with gene_id/gene_name or Name=/product= packed into a
## single "Attributes" text column) -- not a clean gene-by-sample matrix.
## This script parses that format directly rather than pre-cleaning it by
## hand, so the parsing logic itself is visible and checkable.

library(data.table)

raw <- fread("../data/GSE67758_raw_gene_quantifications.csv.gz", sep = "\t", header = TRUE)
cat("Loaded", nrow(raw), "total feature rows (host + pathogen combined, sense + antisense).\n")

# REAL DISCOVERY, DISCLOSED: this file counts sense- and antisense-strand
# reads as SEPARATE rows for the same genomic feature (see the
# "Orientation of counted reads..." column). For a standard gene
# expression count matrix, only the sense-strand counts are the
# biologically meaningful "this gene's own transcript" signal --
# antisense counts reflect antisense transcription, a real but different
# phenomenon, and would double the row count per gene if not filtered out
# here. Confirmed by checking pasilla-style pairs of "sense"/"anti-sense"
# rows sharing identical coordinates in the raw file before filtering.
orientation_col <- "Orientation of counted reads relative to the strand location of the annotation"
cat("Orientation values present:", paste(unique(raw[[orientation_col]]), collapse = ", "), "\n")
raw <- raw[raw[[orientation_col]] == "sense"]
cat("After keeping sense-strand rows only:", nrow(raw), "feature rows.\n")

# Host = human chromosomes (chr1, chr2, ..., chrX, chrY, chrM, plus
# unplaced scaffolds like GL000220.1). Pathogen = the Salmonella SL1344
# chromosome (NC_016810.1) and its three plasmids (NC_017718/019/720.1) --
# listed explicitly, the same real, verified contig list Script 6 of the
# HPC pipeline uses, not inferred.
pathogen_contigs <- c("NC_016810.1", "NC_017718.1", "NC_017719.1", "NC_017720.1")
is_pathogen <- raw$`Sequence name` %in% pathogen_contigs
is_host <- grepl("^chr", raw$`Sequence name`)

cat("Host feature rows:", sum(is_host), "| Pathogen feature rows:", sum(is_pathogen), "\n")

# Extract a usable gene ID/name from the packed Attributes column. THREE
# different real conventions coexist in this one file, verified by direct
# inspection, not assumed to be uniform:
#   - Host rows (GENCODE-style):        "gene_id=...;...;gene_name=...;"
#   - Pathogen CDS rows (RefSeq GFF):    "ID=cds0;Name=YP_...;...;product=...;"
#   - Pathogen sRNA/ncRNA rows:          "locus_tag=ncRNA0164;name=pinT;status=confirmed"
#     (lowercase "name=", "locus_tag=" instead of "ID=" -- this is exactly
#     the row format the paper's headline finding, the sRNA "pinT", uses;
#     an extraction pattern that only looked for "ID=" and "Name=" would
#     silently fail to identify it, so both cases are handled explicitly.)
extract_attr_vectorized <- function(attributes, keys) {
  result <- rep(NA_character_, length(attributes))
  for (key in keys) {
    m <- regmatches(attributes, regexpr(paste0("(^|;)", key, "=[^;]+"), attributes))
    hit <- nchar(m) > 0 & is.na(result)
    result[hit] <- sub(paste0("^;?", key, "="), "", m[hit])
  }
  result
}

host_rows <- raw[is_host]
host_rows[, gene_id := extract_attr_vectorized(Attributes, "gene_id")]
host_rows[, gene_name := extract_attr_vectorized(Attributes, "gene_name")]

pathogen_rows <- raw[is_pathogen]
pathogen_rows[, gene_id := extract_attr_vectorized(Attributes, c("locus_tag", "ID"))]
pathogen_rows[, gene_name := extract_attr_vectorized(Attributes, c("name", "Name"))]
pathogen_rows[, product := extract_attr_vectorized(Attributes, "product")]

pint_check <- pathogen_rows[gene_name == "pinT"]
cat("\npinT (the paper's headline sRNA finding) located:", nrow(pint_check), "row(s), gene_id =",
    if (nrow(pint_check) > 0) pint_check$gene_id[1] else "NOT FOUND -- investigate before trusting the pathogen extraction", "\n")

# This project's 8 real sample columns (see ../data/metadata_8samples.csv,
# `published_column`).
sample_cols <- metadata$published_column
stopifnot(all(sample_cols %in% colnames(raw)))

host_counts <- as.data.frame(host_rows[, c("gene_id", "gene_name", sample_cols), with = FALSE])
host_counts <- host_counts[!duplicated(host_counts$gene_id), ]
rownames(host_counts) <- host_counts$gene_id

pathogen_counts <- as.data.frame(pathogen_rows[, c("gene_id", "gene_name", "product", sample_cols), with = FALSE])
pathogen_counts <- pathogen_counts[!duplicated(pathogen_counts$gene_id), ]
rownames(pathogen_counts) <- pathogen_counts$gene_id

cat("\nHost count matrix:", nrow(host_counts), "genes x", length(sample_cols), "samples\n")
cat("Pathogen count matrix:", nrow(pathogen_counts), "features x", length(sample_cols), "samples\n")

host_matrix <- round(as.matrix(host_counts[, sample_cols]))
pathogen_matrix <- round(as.matrix(pathogen_counts[, sample_cols]))
storage.mode(host_matrix) <- "integer"
storage.mode(pathogen_matrix) <- "integer"

cat("\nHost library sizes:\n"); print(colSums(host_matrix))
cat("\nPathogen library sizes (raw pathogen-mapped read counts per sample):\n"); print(colSums(pathogen_matrix))

cat("\nOK: host_matrix and pathogen_matrix built, columns match metadata's published_column order.\n")
