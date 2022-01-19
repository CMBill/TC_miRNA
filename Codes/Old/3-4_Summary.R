library(tidyverse)
Path <- getwd()

# 表达数据
# exp.select <- read.csv(paste0(Path, "/Data/Rexp_select.csv"), row.names = 1)
# SampleGroup <- read.csv(paste0(Path, './Data/RSampleGroup.csv'), row.names = 1)
# group <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

# 读取三个包的结果
res.DESeq2 <- as.data.frame(read.csv(paste0(Path, '/Data/res_DESeq2.csv'), row.names = 1))
res.edgeR <- as.data.frame(read.csv(paste0(Path, '/Data/res_edgeR.csv'), row.names = 1))
res.limma <- as.data.frame(read.csv(paste0(Path, '/Data/res_limma.csv'), row.names = 1))

# 查看结果
res.summary <- data.frame(row.names = c('Down', 'Stable', 'Up'),
                      DESeq2 = as.integer(table(res.DESeq2$tag)),
                      edgeR = as.integer(table(res.edgeR$tag)),
                      limma = as.integer(table(res.limma$tag))
)

# 分别获取各自DEGs
degs.DESeq2 <- row.names(res.DESeq2[res.DESeq2$tag != 'Stable', ])
degs.edgeR <- row.names(res.edgeR[res.edgeR$tag != 'Stable', ])
degs.limma <- row.names(res.limma[res.limma$tag != 'Stable', ])

degs.all <- intersect(degs.limma, intersect(degs.DESeq2, degs.edgeR))
degs <- merge(res.DESeq2[degs.all, c('log2FoldChange', 'padj')], res.edgeR[degs.all, c('logFC', 'FDR')], by = 'row.names')
degs <- rename(degs, log2FC.DESeq2 = log2FoldChange, FDR.DESeq2 = padj, log2FC.edgeR = logFC, FDR.edgeR = FDR)
row.names(degs) <- degs$Row.names
degs <- degs[, -1]
degs <- merge(degs, res.limma[degs.all, c('logFC', 'adj.P.Val')], by = 'row.names')
degs <- rename(degs, log2FC.limma = logFC, FDR.limma = adj.P.Val)
row.names(degs) <- degs$Row.names
degs <- degs[, -1]

write.csv(res.summary, paste0(Path, '/Data/res_summary.csv'))
write.csv(degs, paste0(Path, '/Data/degs.csv'))