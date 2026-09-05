## 05_volcano_and_top_genes.R
## Volcano plot from Script 03's real DESeq2 host results.

library(ggplot2)
library(ggrepel)

res_host_df$significant <- with(res_host_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
cat("Host genes at padj < 0.05 & |log2FC| > 1:", sum(res_host_df$significant, na.rm = TRUE), "\n")

top_labels <- head(res_host_df[order(res_host_df$padj), ], 10)

p_volcano <- ggplot(res_host_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.4, size = 1) +
  geom_text_repel(data = top_labels, aes(label = gene_name), size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#C00000")) +
  labs(title = "Volcano plot: Salmonella-infected vs mock HeLa (host transcriptome)",
       x = "log2 fold change (infected / mock)", y = "-log10(p-value)",
       color = "padj<0.05 & |log2FC|>1") +
  theme_minimal(base_size = 13)

ggsave("../results/plot_volcano_host.png", p_volcano, width = 7.5, height = 6, dpi = 150)

top20 <- head(res_host_df[order(res_host_df$padj), c("gene_id", "gene_name", "log2FoldChange", "padj")], 20)
write.csv(top20, "../results/top20_host_DE_genes.csv", row.names = FALSE)
cat("\nTop 20 differentially expressed host genes:\n")
print(top20)

cat("\nOK: plot_volcano_host.png and top20_host_DE_genes.csv written to ../results/\n")
