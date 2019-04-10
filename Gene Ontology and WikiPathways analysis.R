#setwd("D:/CityU_research/Auto/Heatmap")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE") 
BiocManager::install("topGO")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("GSEABase")
BiocManager::install("RCy3")
BiocManager::install("rWikiPathways")
BiocManager::install("org.Sc.sgd.db")


library(DOSE)
library(org.Mm.eg.db)## mouse
library(org.Hs.eg.db)## human
library(topGO)
library(clusterProfiler)
library(pathview)
library(org.Sc.sgd.db) # Saccharomyces 
library(rWikiPathways)


keytypes(org.Mm.eg.db) 

data2 <- read.table("GSE78779Expression.txt", header=TRUE)
data$cell_id <- as.character(data$cell_id)
head(data$cell_id, 2)
data_cell_id <- data$cell_id


top <- read.table("GSE78779_top.csv", header=FALSE)
top$V1 <- top$V1 + 1
top$V1[1]
top_cell_id <- data_cell_id[top$V1]
dim(top_cell_id)
neuron_1 <- top_cell_id[1:1997]
neuron_2 <- top_cell_id[1998:3994]
neuron_3 <- top_cell_id[3995:5991]
neuron_4 <- top_cell_id[5992:7988]
neuron_5 <- top_cell_id[7989:9985]
neuron_6 <- top_cell_id[9986:11982]
neuron_7 <- top_cell_id[11983:13979]
neuron_8 <- top_cell_id[13980:15976]
neuron_9 <- top_cell_id[15977:17973]
neuron_10 <- top_cell_id[17974:19970]
#test1 = bitr(data$cell_id, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
#head(test1,2)

test2 = bitr(neuron_1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")


ego_ALL <- enrichGO(gene = test2$ENTREZID, 
                    OrgDb = org.Mm.eg.db,
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #or one of CC,  BP, and  MF
                    pAdjustMethod = "BH", #adjust methods "holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
                    pvalueCutoff = 0.01, 
                    qvalueCutoff = 0.05,
                    readable = TRUE) #Gene ID transfer to gene Symbol
#head(ego_ALL)
dim(ego_ALL)
dim(ego_ALL[ego_ALL$ONTOLOGY=='BP',])
dim(ego_ALL[ego_ALL$ONTOLOGY=='CC',])
dim(ego_ALL[ego_ALL$ONTOLOGY=='MF',])
#barplot(ego_ALL, showCategory=30,title="EnrichmentGO_BP")

ego_BP <- enrichGO(gene = test2$ENTREZID, 
                   OrgDb = org.Mm.eg.db, 
                   #keytype = 'ENTREZID',
                   ont = "BP",
                   pAdjustMethod = "BH", #"holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                   pvalueCutoff = 0.01, 
                   qvalueCutoff = 0.05,
                   readable = TRUE) #Gene ID transfer to gene Symbol
head(ego_BP, 10)
ego_BP
# ego_BP <- setReadable(ego_BP, OrgDb = org.Mm.eg.db)

pdf("GSE60361_neuron_bar_1.pdf", width =10, height = 6)
barplot(ego_BP, showCategory=30,title="EnrichmentGO_BP")
dev.off()

pdf("GSE60361_neuron_cne_1.pdf", width =11, height = 8)
cnetplot(ego_BP,categorySize="pvalue",circular = TRUE, colorEdge = TRUE)
dev.off()

######################################### wikipathways #################################################################################

gene <- test2$ENTREZID
#wikipathways.Sc.gmt <- rWikiPathways::downloadPathwayArchive(organism="Saccharomyces cerevisiae", format = "gmt")
#wp2gene <- read.gmt(wikipathways.Sc.gmt)
wp2gene <- read.gmt("wikipathways-20180810-gmt-Mus_musculus.gmt")
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME


ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
#ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff = 1, qvalueCutoff = 0.5)
ewp <- setReadable(ewp, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
head(ewp)

write.csv(ewp,"GSE60361_neuron_1.csv")
