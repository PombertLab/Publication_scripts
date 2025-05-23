#!/usr/bin/Rscript
library(ComplexHeatmap)
library(RColorBrewer)
library(methods)
library(circlize)
colors <- colorRamp2(c(0, 20, 40, 60, 80, 100), c("white", "yellow", "lightblue", "blue", "magenta", "red"))
ht_global_opt(heatmap_row_names_gp = gpar(fontsize = 1, fontface = "italic"), heatmap_column_names_gp = gpar(fontsize = 1), heatmap_column_title_gp = gpar(fontsize = 1))
BATHYvsBATHY <- read.csv("BATHYvsBATHY.matrix", header=TRUE)
rownames(BATHYvsBATHY) <- BATHYvsBATHY[,1]
colnames(BATHYvsBATHY)
data_BATHYvsBATHY <- data.matrix(BATHYvsBATHY[,2:ncol(BATHYvsBATHY)])
ht_BATHYvsBATHY = Heatmap(data_BATHYvsBATHY, name = "BATHYvsBATHY", width = unit(57, "mm"), show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lty = 1, lwd = 0.01),column_title = "BATHYvsBATHY", col = colors)
class(ht_BATHYvsBATHY)
BATHYvsCCMP <- read.csv("BATHYvsCCMP.matrix", header=TRUE)
rownames(BATHYvsCCMP) <- BATHYvsCCMP[,1]
colnames(BATHYvsCCMP)
data_BATHYvsCCMP <- data.matrix(BATHYvsCCMP[,2:ncol(BATHYvsCCMP)])
ht_BATHYvsCCMP = Heatmap(data_BATHYvsCCMP, name = "BATHYvsCCMP", width = unit(60, "mm"), show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lty = 1, lwd = 0.01),column_title = "BATHYvsCCMP", col = colors)
class(ht_BATHYvsCCMP)
BATHYvsCOCCO <- read.csv("BATHYvsCOCCO.matrix", header=TRUE)
rownames(BATHYvsCOCCO) <- BATHYvsCOCCO[,1]
colnames(BATHYvsCOCCO)
data_BATHYvsCOCCO <- data.matrix(BATHYvsCOCCO[,2:ncol(BATHYvsCOCCO)])
ht_BATHYvsCOCCO = Heatmap(data_BATHYvsCOCCO, name = "BATHYvsCOCCO", width = unit(87, "mm"), show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lty = 1, lwd = 0.01),column_title = "BATHYvsCOCCO", col = colors)
class(ht_BATHYvsCOCCO)
BATHYvsMCOMMODA <- read.csv("BATHYvsMCOMMODA.matrix", header=TRUE)
rownames(BATHYvsMCOMMODA) <- BATHYvsMCOMMODA[,1]
colnames(BATHYvsMCOMMODA)
data_BATHYvsMCOMMODA <- data.matrix(BATHYvsMCOMMODA[,2:ncol(BATHYvsMCOMMODA)])
ht_BATHYvsMCOMMODA = Heatmap(data_BATHYvsMCOMMODA, name = "BATHYvsMCOMMODA", width = unit(51, "mm"), show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lty = 1, lwd = 0.01),column_title = "BATHYvsMCOMMODA", col = colors)
class(ht_BATHYvsMCOMMODA)
BATHYvsMPUSI <- read.csv("BATHYvsMPUSI.matrix", header=TRUE)
rownames(BATHYvsMPUSI) <- BATHYvsMPUSI[,1]
colnames(BATHYvsMPUSI)
data_BATHYvsMPUSI <- data.matrix(BATHYvsMPUSI[,2:ncol(BATHYvsMPUSI)])
ht_BATHYvsMPUSI = Heatmap(data_BATHYvsMPUSI, name = "BATHYvsMPUSI", width = unit(60, "mm"), show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lty = 1, lwd = 0.01),column_title = "BATHYvsMPUSI", col = colors)
class(ht_BATHYvsMPUSI)
BATHYvsOLUCI <- read.csv("BATHYvsOLUCI.matrix", header=TRUE)
rownames(BATHYvsOLUCI) <- BATHYvsOLUCI[,1]
colnames(BATHYvsOLUCI)
data_BATHYvsOLUCI <- data.matrix(BATHYvsOLUCI[,2:ncol(BATHYvsOLUCI)])
ht_BATHYvsOLUCI = Heatmap(data_BATHYvsOLUCI, name = "BATHYvsOLUCI", width = unit(63, "mm"), show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lty = 1, lwd = 0.01),column_title = "BATHYvsOLUCI", col = colors)
class(ht_BATHYvsOLUCI)
BATHYvsOTAURI <- read.csv("BATHYvsOTAURI.matrix", header=TRUE)
rownames(BATHYvsOTAURI) <- BATHYvsOTAURI[,1]
colnames(BATHYvsOTAURI)
data_BATHYvsOTAURI <- data.matrix(BATHYvsOTAURI[,2:ncol(BATHYvsOTAURI)])
ht_BATHYvsOTAURI = Heatmap(data_BATHYvsOTAURI, name = "BATHYvsOTAURI", width = unit(66, "mm"), show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lty = 1, lwd = 0.01),column_title = "BATHYvsOTAURI", col = colors)
class(ht_BATHYvsOTAURI)
pdf(file="BATHY.pdf", useDingbats=FALSE, width=23, height=5)
ht_list = ht_BATHYvsBATHY + ht_BATHYvsCCMP + ht_BATHYvsCOCCO + ht_BATHYvsMCOMMODA + ht_BATHYvsMPUSI + ht_BATHYvsOLUCI + ht_BATHYvsOTAURI 
draw(ht_list, gap = unit(0.75, "mm"), heatmap_legend_side = "bottom")
dev.off()
