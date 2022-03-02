
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
genes.table <- as.data.frame(table(miRNA.targets$genesymbol))
genes <- as.character(genes.table$Var1)

