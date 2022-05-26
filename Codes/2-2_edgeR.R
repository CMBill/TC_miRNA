# edgeR
library(edgeR)
Path <- getwd()

s.count <- read.csv(paste0(Path, '/Data/s_count.csv'), row.names = 1, check.names = FALSE)
SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'), row.names = 1)
group <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

designMat <- model.matrix(~0+group)
row.names(designMat) <- row.names(SampleGroup)
dgel <- DGEList(counts = s.count, group = group)
dgel.norm <- calcNormFactors(dgel, method = 'TMM')
dgel.est <- estimateGLMCommonDisp(dgel.norm, designMat)
dgel.est <- estimateGLMTrendedDisp(dgel.est, designMat)
dgel.est <- estimateGLMTagwiseDisp(dgel.est, designMat)
fit <- glmFit(dgel.est, designMat)
lrt <- glmLRT(fit, contrast = c(-1, 1))
# 筛选
res.edgeR <- as.data.frame(topTags(lrt, n = nrow(dgel$counts)))

res.edgeR$tag <- ifelse(res.edgeR$FDR < 0.05 & abs(res.edgeR$logFC) >= 1, 
                        ifelse(res.edgeR$logFC> 1,'Up', 'Down'), 
                        'Stable')
write.csv(res.edgeR, paste0(Path, '/Data/res_edgeR_rpm.csv'))
