## 09_compare_hpc_pipeline_to_published.R
## Requires pipeline_raw_reads_to_counts_output/counts/ (Script 7's real
## output, copied down from the HPC cluster -- see README's "Downstream
## handoff" section).
##
## Compares this project's own STAR/featureCounts-derived host counts
## against the published DESeq2 results from Script 03.
##
## REAL, VERIFIED ID-FORMAT MISMATCH, DISCLOSED (not silently patched):
## the published matrix's gene_id column carries a GENCODE-style version
## suffix (e.g. "ENSG00000166145.10"), while this pipeline's own Ensembl
## GTF (verified directly by inspecting its first bytes) omits the version
## from the gene_id attribute entirely ("ENSG00000279928", with the "10"
## stored separately as "gene_version"). Stripping the suffix from the
## published IDs before joining is the correct fix here, not a Version
## mismatch worth pretending doesn't exist.

host_path <- "../pipeline_raw_reads_to_counts_output/counts/cohort_host_counts.txt"

if (!file.exists(host_path)) {
  stop("Run Script 7 on the HPC cluster first, then copy counts/cohort_host_counts.txt into ",
       "../pipeline_raw_reads_to_counts_output/counts/ -- see README's Downstream handoff section.")
}

hpc_host <- read.table(host_path, header = TRUE, skip = 1, row.names = "Geneid", check.names = FALSE)
hpc_host <- hpc_host[, 6:ncol(hpc_host), drop = FALSE]
colnames(hpc_host) <- sub(".*/(SRR[0-9]+)_.*", "\\1", colnames(hpc_host))

cat("HPC pipeline's own host counts:", nrow(hpc_host), "genes x", ncol(hpc_host), "samples\n")

srr_to_condition <- setNames(metadata$condition, metadata$srr)
hpc_condition <- srr_to_condition[colnames(hpc_host)]
stopifnot(!any(is.na(hpc_condition)))

cpm <- sweep(hpc_host, 2, colSums(hpc_host) / 1e6, "/")
hpc_infected_mean <- rowMeans(cpm[, hpc_condition == "infected", drop = FALSE])
hpc_mock_mean <- rowMeans(cpm[, hpc_condition == "mock", drop = FALSE])
hpc_log2fc <- log2((hpc_infected_mean + 1) / (hpc_mock_mean + 1))

# Strip the GENCODE-style version suffix from the published gene IDs so
# they match this pipeline's own unversioned Ensembl gene_id values.
published_ids_stripped <- sub("\\.[0-9]+$", "", res_host_df$gene_id)
published_log2fc <- setNames(res_host_df$log2FoldChange, published_ids_stripped)

common_genes <- intersect(names(hpc_log2fc), names(published_log2fc))
cat("Genes in common between HPC pipeline output and published matrix (after version-stripping):",
    length(common_genes), "\n")

comparison <- data.frame(
  gene_id = common_genes,
  published_log2FC = published_log2fc[common_genes],
  hpc_log2FC = hpc_log2fc[common_genes]
)
comparison <- comparison[stats::complete.cases(comparison), ]

corr <- cor(comparison$published_log2FC, comparison$hpc_log2FC, method = "spearman")
cat(sprintf("\nSpearman correlation (published DESeq2 log2FC vs HPC pipeline's own log2FC): %.3f\n", corr))
cat("(Modest cohort -- 8 samples, single representative HPC run per sample -- expect a real but\n")
cat(" imperfect correlation, same caveat as this curriculum's other gold-standard comparisons.)\n")

library(ggplot2)
p_compare <- ggplot(comparison, aes(published_log2FC, hpc_log2FC)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "#C00000") +
  labs(title = "Published DESeq2 log2FC vs this project's own HPC pipeline (host side)",
       subtitle = sprintf("Spearman rho = %.3f, n = %d genes", corr, nrow(comparison)),
       x = "Published (GSE67758 supplementary matrix) log2FC",
       y = "This project's own STAR/featureCounts log2FC") +
  theme_minimal(base_size = 13)

dir.create("../results", showWarnings = FALSE)
ggsave("../results/plot_hpc_vs_published_comparison.png", p_compare, width = 7, height = 6, dpi = 150)
write.csv(comparison, "../results/hpc_vs_published_comparison.csv", row.names = FALSE)

cat("\nOK: plot_hpc_vs_published_comparison.png and hpc_vs_published_comparison.csv written to ../results/\n")
