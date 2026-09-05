## 03_differential_expression_deseq2.R
## Host-side DESeq2: infected vs mock, controlling for the real 4h/24h
## timepoint covariate (design: ~ timepoint_h + condition). Real raw
## integer counts from the published matrix -- valid DESeq2 input directly.

library(DESeq2)

col_data <- metadata[colnames(host_matrix), ]
stopifnot(identical(rownames(col_data), colnames(host_matrix)))

dds_host <- DESeqDataSetFromMatrix(countData = host_matrix,
                                    colData = col_data,
                                    design = ~ timepoint_h + condition)

keep <- rowSums(counts(dds_host)) >= 10
cat("Filtering:", sum(!keep), "of", length(keep), "host genes dropped (rowSum < 10 across all 8 samples)\n")
dds_host <- dds_host[keep, ]

dds_host <- DESeq(dds_host)
print(resultsNames(dds_host))
res_host <- results(dds_host, name = grep("^condition_", resultsNames(dds_host), value = TRUE))

res_host_df <- as.data.frame(res_host)
res_host_df$gene_id <- rownames(res_host_df)
res_host_df$gene_name <- host_counts[res_host_df$gene_id, "gene_name"]
res_host_df <- res_host_df[order(res_host_df$padj), ]

cat("\nDESeq2 summary (Salmonella-infected vs mock HeLa, controlling for timepoint):\n")
summary(res_host)

cat("\nTop 15 host genes by adjusted p-value:\n")
print(head(res_host_df[, c("gene_id", "gene_name", "log2FoldChange", "pvalue", "padj")], 15))

dir.create("../results", showWarnings = FALSE)
write.csv(res_host_df, "../results/deseq2_host_infected_vs_mock.csv", row.names = FALSE)
saveRDS(dds_host, "../results/dds_host.rds")
