# DESeq2
library(DESeq2)
Path <- getwd()

# counts
count.select <- read.csv(paste0(Path, '/Data/selected_count.csv'), row.names = 1)
SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'), row.names = 1)

SampleGroup[, 'Group'] <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))
row.names(SampleGroup) <- make.names(row.names(SampleGroup))

dds <- DESeqDataSetFromMatrix(count.select, colData=SampleGroup, design=~Group)
des.dds <- DESeq(dds)
res.dds <- results(des.dds)
res.dds <- as.data.frame(res.dds)

res.dds$tag <- ifelse(res.dds$padj < 0.05 & abs(res.dds$log2FoldChange) >= 1,
                      ifelse(res.dds$log2FoldChange> 1, 'Up', 'Down'),
                      'Stable')

write.csv(res.dds, paste0(Path, '/Data/res_DESeq2.csv'))

# rpm
rpm.select <- read.csv(paste0(Path, '/Data/selected_rpm.csv'), row.names = 1)

dds.rpm <- DESeqDataSetFromMatrix(round(rpm.select), colData=SampleGroup, design=~Group)
des.dds.rpm <- DESeq(dds.rpm)
res.dds.rpm <- results(des.dds.rpm)
res.dds.rpm <- as.data.frame(res.dds.rpm)

res.dds.rpm$tag <- ifelse(res.dds.rpm$padj < 0.05 & abs(res.dds.rpm$log2FoldChange) >= 1, 
                          ifelse(res.dds.rpm$log2FoldChange> 1, 'Up', 'Down'), 
                          'Stable')

write.csv(res.dds.rpm, paste0(Path, '/Data/res_DESeq2_rpm.csv'))
