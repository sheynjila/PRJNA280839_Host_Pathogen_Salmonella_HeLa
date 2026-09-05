## 04_exploratory_overview.R
## PCA and sample-correlation heatmap on the host side (VST-transformed).

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

vsd_host <- vst(dds_host, blind = TRUE)

pca_data <- plotPCA(vsd_host, intgroup = c("condition", "timepoint_h"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = timepoint_h)) +
  geom_point(size = 4, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c(mock = "#2E74B5", infected = "#C00000")) +
  labs(title = "PCA: Salmonella-infected vs mock HeLa cells",
       subtitle = "Shape shows timepoint (4h / 24h post-infection)",
       color = "Condition", shape = "Timepoint") +
  theme_minimal(base_size = 13)

ggsave("../results/plot_pca_host.png", p_pca, width = 7, height = 5.5, dpi = 150)

sample_dist <- dist(t(assay(vsd_host)))
sample_dist_matrix <- as.matrix(sample_dist)

png("../results/plot_sample_correlation_heatmap.png", width = 900, height = 800, res = 150)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distance (host, VST-transformed counts)")
dev.off()

cat("OK: plot_pca_host.png and plot_sample_correlation_heatmap.png written to ../results/\n")
cat("Expect PC1 or PC2 to separate by TIMEPOINT (4h vs 24h) as strongly as by condition --\n")
cat("infection response is a time-dependent process, not a single fixed state.\n")
