#!/usr/bin/Rscript
library(ComplexHeatmap)
library(RColorBrewer)
library(methods)
message ("Plotting COCCOvsOLUCI.matrix...")
colors <- colorRampPalette(c("white", "blue", "magenta"))(n = 300)
ht_global_opt(heatmap_row_names_gp = gpar(fontsize = 8, fontface = "italic"), heatmap_column_names_gp = gpar(fontsize = 8), heatmap_column_title_gp = gpar(fontsize = 12))
pdf(file="COCCOvsOLUCI.matrix.pdf", useDingbats=FALSE, width=5, height=4)
COCCOvsOLUCI <- read.csv("COCCOvsOLUCI.matrix", header=TRUE)
rownames(COCCOvsOLUCI) <- COCCOvsOLUCI[,1]
colnames(COCCOvsOLUCI)
data_COCCOvsOLUCI <- data.matrix(COCCOvsOLUCI[,2:ncol(COCCOvsOLUCI)])
ht_COCCOvsOLUCI = Heatmap(data_COCCOvsOLUCI, name = "COCCOvsOLUCI", width = unit(63, "mm"), cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lty = 1, lwd = 1), column_title = "COCCOvsOLUCImatrix", col = colors)
class(ht_COCCOvsOLUCI)
draw(ht_COCCOvsOLUCI, heatmap_legend_side = "right")
dev.off()
