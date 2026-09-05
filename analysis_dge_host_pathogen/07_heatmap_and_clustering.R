## 07_heatmap_and_clustering.R
## Top 40 DE host genes, row Z-scores, Ward's-method clustering.

library(pheatmap)
library(DESeq2)
library(RColorBrewer)

top40_ids <- head(res_host_df[order(res_host_df$padj), "gene_id"], 40)
vst_mat <- assay(vst(dds_host, blind = TRUE))[top40_ids, ]
rownames(vst_mat) <- host_counts[rownames(vst_mat), "gene_name"]

annotation_col <- data.frame(
  condition = metadata[colnames(vst_mat), "condition"],
  timepoint_h = metadata[colnames(vst_mat), "timepoint_h"],
  row.names = colnames(vst_mat)
)

png("../results/plot_heatmap_top40_host.png", width = 1000, height = 1100, res = 150)
pheatmap(vst_mat,
         scale = "row",
         clustering_method = "ward.D2",
         annotation_col = annotation_col,
         annotation_colors = list(
           condition = c(mock = "#2E74B5", infected = "#C00000"),
           timepoint_h = c(`4` = "#B08D57", `24` = "#7A5200")
         ),
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
         main = "Top 40 DE host genes (Salmonella-infected vs mock)",
         show_rownames = TRUE, fontsize_row = 6)
dev.off()

cat("OK: plot_heatmap_top40_host.png written to ../results/\n")
