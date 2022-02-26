# DESeq2
library(DESeq2)
Path <- getwd()

s.count <- read.csv(paste0(Path, '/Data/s_count.csv'), row.names = 1)
SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'), row.names = 1)
SampleGroup[, 'Group'] <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))
row.names(SampleGroup) <- make.names(row.names(SampleGroup))

dds <- DESeqDataSetFromMatrix(round(s.count), colData=SampleGroup, design=~Group)
des.dds <- DESeq(dds)
res.DESeq2 <- results(des.dds)
res.DESeq2 <- as.data.frame(res.DESeq2)

res.DESeq2$tag <- ifelse(res.DESeq2$padj < 0.05 & abs(res.DESeq2$log2FoldChange) >= 1, 
                         ifelse(res.DESeq2$log2FoldChange> 1, 'Up', 'Down'), 
                         'Stable')

write.csv(res.DESeq2, paste0(Path, '/Data/res_DESeq2_rpm.csv'))
