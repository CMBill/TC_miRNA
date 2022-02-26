library(DESeq2)
Path <- getwd()

sub_count <- read.csv(paste0(Path, '/Data/sub_count.csv'), row.names = 1)
sub_rpm <- read.csv(paste0(Path, '/Data/sub_rpm.csv'), row.names = 1)
sub_group <- read.csv(paste0(Path, './Data/sub_group.csv'))

s_ptc <- sub_group[sub_group$Group == 'PTC', ]
s_fvptc <- sub_group[sub_group$Group == 'FV-PTC', ]
s_nor <- sub_group[sub_group$Group == 'normal', ]

s.group <- rbind(s_ptc, s_fvptc)
row.names(s.group) <- make.names(s.group$EntitiesId)
s.group <- s.group[, -1, drop = FALSE]
s_group <- data.frame()

s.count <- sub_count[, colnames(sub_count) %in% row.names(s.group)]
s_group <- s.group[colnames(s.count), , drop = FALSE]

s_group$Group <- factor(s_group$Group, levels = c('PTC', 'FV-PTC'), labels = c('PTC', 'FVPTC'))

dds <- DESeqDataSetFromMatrix(round(s.count), colData=s_group, design=~Group)
des.dds <- DESeq(dds)
res.DESeq2 <- results(des.dds)
res.DESeq2 <- as.data.frame(res.DESeq2)

res.DESeq2$tag <- ifelse(res.DESeq2$padj < 0.05 & abs(res.DESeq2$log2FoldChange) >= 1, 
                         ifelse(res.DESeq2$log2FoldChange> 1, 'Up', 'Down'), 
                         'Stable')

#edgeR
group <- factor(s_group$Group, levels = c('PTC', 'FV-PTC'), labels = c('PTC', 'FVPTC'))
designMat <- model.matrix(~0+group)
row.names(designMat) <- row.names(s_group)

# 标准化
dgel <- DGEList(counts = s.count, group = group)
dgel.norm <- calcNormFactors(dgel, method = 'TMM')

# Disperison估算
dgel.est <- estimateGLMCommonDisp(dgel.norm, designMat)
dgel.est <- estimateGLMTrendedDisp(dgel.est, designMat)
dgel.est <- estimateGLMTagwiseDisp(dgel.est, designMat)

# 拟合
fit <- glmFit(dgel.est, designMat)
lrt <- glmLRT(fit, contrast = c(-1, 1))

# 筛选
res.edgeR <- as.data.frame(topTags(lrt, n = nrow(dgel$counts)))

res.edgeR$tag <- ifelse(res.edgeR$FDR < 0.05 & abs(res.edgeR$logFC) >= 1, 
                        ifelse(res.edgeR$logFC> 1,'Up', 'Down'), 
                        'Stable')
