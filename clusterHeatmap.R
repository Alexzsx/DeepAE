
library("ComplexHeatmap")
library(circlize)
setwd("D:/CityU Researches/Auto/Heatmap/")



png(file="OriginalGSE60361_3.png")
Original_data = read.table("OrigianlGSE60361.csv", header=F,sep = ",")

mat = log1p(Original_data[])
#column_title = "Original data"
Heatmap(mat, column_title = "Samples", name = "",
        row_dend_reorder = FALSE, show_column_names = FALSE,
        row_title = "Genes", column_title_side = "bottom")
dev.off()


png(file="DeepAEGSE60361_3.png")
Original_data = read.table("DeepAEGSE60361.csv", header=F,sep = ",")

mat = log1p(Original_data[])
#column_title = "DeepAE"
Heatmap(mat, column_title = "Samples", name = "", row_dend_reorder = FALSE, 
        show_column_names = FALSE, 
        row_title = "Genes", column_title_side = "bottom")
dev.off()

png(file="GSE60361_svd2_3.png")
Original_data = read.table("GSE60361_svd2.csv", header=F,sep = ",")
Original_data[is.na(Original_data)]=0.001
mat = log1p(Original_data[]*10)
#column_title = "SVD"
Heatmap(mat, column_title = "Samples", name = "", row_dend_reorder = FALSE, 
        show_column_names = FALSE, 
        row_title = "Genes", column_title_side = "bottom")
dev.off()

png(file="GSE60361_ksvd2_3.png")
Original_data = read.table("GSE60361_ksvd2.csv", header=F,sep = ",")
Original_data[is.na(Original_data)]=0.001
mat = log1p(Original_data[]*10)
#column_title = "k-SVD"
Heatmap(mat, column_title = "Samples", name = "", row_dend_reorder = FALSE, 
        show_column_names = FALSE, 
        row_title = "Genes", column_title_side = "bottom")
dev.off()

png(file="GSE60361_snmf2_3.png")
Original_data = read.table("GSE60361_snmf2.csv", header=F,sep = ",")
#column_title = "sNMF"
Original_data[is.na(Original_data)]=0.001
mat = log1p(Original_data[]*10)

Heatmap(mat, column_title = "Samples", name = "", row_dend_reorder = FALSE, 
        show_column_names = FALSE, 
        row_title = "Genes", column_title_side = "bottom")
dev.off()


png(file="GSE60361_smaf2_3.png")
Original_data = read.table("GSE60361_smaf2.csv", header=F,sep = ",")
Original_data[is.na(Original_data)]=0.001
mat = log1p(Original_data[]*10)
#column_title = "CS-SMAF"
Heatmap(mat, column_title = "Samples", name = "", row_dend_reorder = FALSE, 
        show_column_names = FALSE, 
        row_title = "Genes", column_title_side = "bottom")
dev.off()






