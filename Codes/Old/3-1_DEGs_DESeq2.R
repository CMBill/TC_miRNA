# DESeq2
library(DESeq2)
Path <- getwd()

# 不需要标准化，直接导入
exp_select <- read.csv(paste0(Path, '/Data/Rexp_select.csv'), row.names = 1)
SampleGroup <- read.csv(paste0(Path, './Data/RSampleGroup.csv'), row.names = 1)
SampleGroup[, 'Group'] <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

dds <- DESeqDataSetFromMatrix(exp_select, colData=SampleGroup, design=~Group)
des.dds <- DESeq(dds)
res.dds <- results(des.dds)
res.dds <- as.data.frame(res.dds)

res.dds$tag <- ifelse(res.dds$padj < 0.05 & abs(res.dds$log2FoldChange) >= 1,
                      ifelse(res.dds$log2FoldChange> 1, 'Up', 'Down'),
                      'Stable')

write.csv(res.dds, paste0(Path, '/Data/res_DESeq2.csv'))
