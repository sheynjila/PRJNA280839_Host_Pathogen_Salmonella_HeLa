## 08_publication_dashboard.R
## One composite 4-panel figure (PCA, volcano, pinT showcase, p-value
## histogram) via patchwork.

library(ggplot2)
library(patchwork)
library(DESeq2)

vsd_host <- vst(dds_host, blind = TRUE)
pca_data <- plotPCA(vsd_host, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p1 <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(mock = "#2E74B5", infected = "#C00000")) +
  labs(title = "A. PCA (host)", x = paste0("PC1 (", percentVar[1], "%)"), y = paste0("PC2 (", percentVar[2], "%)")) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

p2 <- ggplot(res_host_df, aes(log2FoldChange, -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#C00000")) +
  labs(title = "B. Volcano (host)") + theme_minimal(base_size = 11) + theme(legend.position = "none")

p3 <- ggplot(pint_expr, aes(interaction(condition, timepoint_h, sep = "@"), raw_count, color = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(mock = "#2E74B5", infected = "#C00000")) +
  labs(title = "C. pinT sRNA (pathogen)", x = NULL, y = "Raw counts") +
  theme_minimal(base_size = 9) + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p4 <- ggplot(res_host_df, aes(pvalue)) +
  geom_histogram(bins = 40, fill = "#2E74B5", color = "white") +
  labs(title = "D. Host p-value distribution", x = "p-value", y = "Gene count") +
  theme_minimal(base_size = 11)

dashboard <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "PRJNA280839: Salmonella-HeLa dual RNA-seq (Westermann et al. 2016)",
                   theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave("../results/plot_publication_dashboard.png", dashboard, width = 11, height = 9, dpi = 150)
cat("OK: plot_publication_dashboard.png written to ../results/\n")
