library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)

Path <- getwd()

res.DESeq2 <- read.csv(paste0(Path, '/Data/res_DESeq2_rpm.csv'), row.names = 1)
res.edgeR <- read.csv(paste0(Path, '/Data/res_edgeR_rpm.csv'), row.names = 1)
res.sum <- data.frame(row.names = c('Down', 'Stable', 'Up'),
                      DESeq2 = as.integer(table(res.DESeq2$tag)),
                      edgeR = as.integer(table(res.edgeR$tag))
                      )

degs.DESeq2 <- row.names(res.DESeq2[res.DESeq2$tag != 'Stable', ])
degs.DESeq2.up <- row.names(res.DESeq2[res.DESeq2$tag == 'Up', ])
degs.DESeq2.down <- row.names(res.DESeq2[res.DESeq2$tag == 'Down', ])

degs.edgeR <- row.names(res.edgeR[res.edgeR$tag != 'Stable', ])
degs.edgeR.up <- row.names(res.edgeR[res.edgeR$tag == 'Up', ])
degs.edgeR.down <- row.names(res.edgeR[res.edgeR$tag == 'Down', ])

degs <- union(degs.DESeq2, degs.edgeR)
degs.up <- union(degs.DESeq2.up, degs.edgeR.up)
degs.down <- union(degs.DESeq2.down, degs.edgeR.down)

# 读取靶文件
miRNA.targets <- read.csv('./Data/miRWalk_miRNA_Targets.csv')
miRNA.targets1 <- read.csv('./Data/miRWalk_miRNA_Targets1.csv')

miRNA.gene1 <- miRNA.targets1[, c('锘縨irnaid', 'genesymbol')]

genes.table <- as.data.frame(table(miRNA.targets$genesymbol))
genes <- as.character(genes.table$Var1)
DEGs_entrez_id <- mapIds(x = org.Hs.eg.db,column = "ENTREZID", keys = genes, keytype = "SYMBOL")
enrich_go_BP <- enrichGO(gene = DEGs_entrez_id,OrgDb = "org.Hs.eg.db",ont = "BP",pvalueCutoff = 0.05)
enrich_go_MF <- enrichGO(gene = DEGs_entrez_id,OrgDb = "org.Hs.eg.db",ont = "MF",pvalueCutoff = 0.05)
enrich_go_CC <- enrichGO(gene = DEGs_entrez_id,OrgDb = "org.Hs.eg.db",ont = "CC",pvalueCutoff = 0.05)
enrich_go_KEGG <- enrichKEGG(gene = DEGs_entrez_id,organism = "hsa",keyType = "kegg",pvalueCutoff = 0.05)

BP <- dotplot(enrich_go_BP,showCategory=20)
MF <- dotplot(enrich_go_MF,showCategory=20)
CC <- dotplot(enrich_go_CC,showCategory=20)
KEGG <- dotplot(enrich_go_KEGG)

miRNA.table1 <- as.data.frame(table(miRNA.targets1$锘縨irnaid))
miRNA.table1.order <- miRNA.table1[order(miRNA.table1$Freq, decreasing = T),]
miRNA.top7 <- miRNA.table1.order[1:7, ]
miRNA.t7 <- as.character(miRNA.top7$Var1)

genes.table1 <- as.data.frame(table(miRNA.targets1$genesymbol))
genes1 <- as.character(genes.table1$Var1)

DEGs_entrez_id1 <- mapIds(x = org.Hs.eg.db,column = "ENTREZID", keys = genes1, keytype = "SYMBOL")
enrich_go_BP1 <- enrichGO(gene = DEGs_entrez_id1,OrgDb = "org.Hs.eg.db",ont = "BP",pvalueCutoff = 0.05)
enrich_go_MF1 <- enrichGO(gene = DEGs_entrez_id1,OrgDb = "org.Hs.eg.db",ont = "MF",pvalueCutoff = 0.05)
enrich_go_CC1 <- enrichGO(gene = DEGs_entrez_id1,OrgDb = "org.Hs.eg.db",ont = "CC",pvalueCutoff = 0.05)
enrich_go_KEGG1 <- enrichKEGG(gene = DEGs_entrez_id1,organism = "hsa",keyType = "kegg",pvalueCutoff = 0.05)

BP1 <- dotplot(enrich_go_BP1,showCategory=10)
MF1 <- dotplot(enrich_go_MF1,showCategory=10)
CC1 <- dotplot(enrich_go_CC1,showCategory=10)
KEGG1 <- dotplot(enrich_go_KEGG1)

miRNA.targets1.top7 <- miRNA.gene1[miRNA.gene1$锘縨irnaid %in% miRNA.t7,]
miRNA.gene1.table.t7 <- as.data.frame(table(miRNA.targets1.top7$genesymbol))
genes1.t7 <- as.character(miRNA.gene1.table.t7$Var1)

DEGs_entrez_id1.t7 <- mapIds(x = org.Hs.eg.db,column = "ENTREZID", keys = genes1.t7, keytype = "SYMBOL")
enrich_go_BP1.t7 <- enrichGO(gene = DEGs_entrez_id1,OrgDb = "org.Hs.eg.db",ont = "BP",pvalueCutoff = 0.05)
enrich_go_MF1.t7 <- enrichGO(gene = DEGs_entrez_id1,OrgDb = "org.Hs.eg.db",ont = "MF",pvalueCutoff = 0.05)
enrich_go_CC1.t7 <- enrichGO(gene = DEGs_entrez_id1,OrgDb = "org.Hs.eg.db",ont = "CC",pvalueCutoff = 0.05)
enrich_go_KEGG1.t7 <- enrichKEGG(gene = DEGs_entrez_id1,organism = "hsa",keyType = "kegg",pvalueCutoff = 0.05)

BP1.t7 <- dotplot(enrich_go_BP1.t7,showCategory=10)
MF1.t7 <- dotplot(enrich_go_MF1.t7,showCategory=10)
CC1.t7 <- dotplot(enrich_go_CC1.t7,showCategory=10)
KEGG1.t7 <- dotplot(enrich_go_KEGG1.t7)
