#!/usr/bin/Rscript
library(ComplexHeatmap)
library(RColorBrewer)
library(methods)
message ("Plotting BATHYvsCCMP.matrix...")
colors <- colorRampPalette(c("white", "blue", "magenta"))(n = 300)
ht_global_opt(heatmap_row_names_gp = gpar(fontsize = 8, fontface = "italic"), heatmap_column_names_gp = gpar(fontsize = 8), heatmap_column_title_gp = gpar(fontsize = 12))
pdf(file="BATHYvsCCMP.matrix.pdf", useDingbats=FALSE, width=5, height=4)
BATHYvsCCMP <- read.csv("BATHYvsCCMP.matrix", header=TRUE)
rownames(BATHYvsCCMP) <- BATHYvsCCMP[,1]
colnames(BATHYvsCCMP)
data_BATHYvsCCMP <- data.matrix(BATHYvsCCMP[,2:ncol(BATHYvsCCMP)])
ht_BATHYvsCCMP = Heatmap(data_BATHYvsCCMP, name = "BATHYvsCCMP", width = unit(60, "mm"), cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lty = 1, lwd = 1), column_title = "BATHYvsCCMPmatrix", col = colors)
class(ht_BATHYvsCCMP)
draw(ht_BATHYvsCCMP, heatmap_legend_side = "right")
dev.off()
