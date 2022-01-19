# edgeR
library(edgeR)
Path <- getwd()

# 导入数据
exp_select <- read.csv(paste0(Path, '/Data/Rexp_select.csv'), row.names = 1)
SampleGroup <- read.csv(paste0(Path, './Data/RSampleGroup.csv'), row.names = 1)
group <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

# 标准化
dgel <- DGEList(counts = exp_select, group = group)
dgel_norm <- calcNormFactors(dgel, method = 'TMM')

# 设计分组信息
designMat <- model.matrix(~0+group)
row.names(designMat) <- row.names(SampleGroup)

# Disperison估算
dgel.est <- estimateGLMCommonDisp(dgel_norm, designMat)
dgel.est <- estimateGLMTrendedDisp(dgel.est, designMat)
dgel.est <- estimateGLMTagwiseDisp(dgel.est, designMat)

# 拟合
fit <- glmFit(dgel.est, designMat)
lrt <- glmLRT(fit, contrast = c(-1, 1))

# 筛选
res.dgel_norm <- as.data.frame(topTags(lrt, n = nrow(dgel$counts)))

res.dgel_norm$tag <- ifelse(res.dgel_norm$FDR < 0.05 & abs(res.dgel_norm$logFC) >= 1,
                            ifelse(res.dgel_norm$logFC> 1, 'Up', 'Down'),
                            'Stable')

write.csv(res.dgel_norm, paste0(Path, '/Data/res_edgeR.csv'))