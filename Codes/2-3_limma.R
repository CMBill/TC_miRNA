# limma
library(limma)
Path <- getwd()

# 导入数据
s.count <- read.csv(paste0(Path, '/Data/s_count.csv'), row.names = 1)
SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'), row.names = 1)
group <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

# 设计分组信息
designMat <- model.matrix(~0+group)
row.names(designMat) <- row.names(SampleGroup)

# 标准化
dgel <- DGEList(counts = s.count, group = group)
dgel.norm <- calcNormFactors(dgel, method = 'TMM')

dgel.v <- voom(dgel.norm, designMat, plot = TRUE, normalize = 'quantile')
fit.lm <- lmFit(dgel.v, designMat)

# 定义比较分组关系
contrasts.lm <- paste(rev(colnames(designMat)),collapse = "-")
cont.matrix <- makeContrasts(contrasts = contrasts.lm, levels = designMat)

# 比较每个基因
fit.lm2 <- contrasts.fit(fit.lm, cont.matrix)
fit.lm2 <- eBayes(fit.lm2)

# 获取结论
fit.lm2.res <- na.omit(topTable(fit.lm2, coef = contrasts.lm, n = Inf))

# 筛选
res.lm <- as.data.frame(fit.lm2.res)

res.lm$tag <- ifelse(res.lm$adj.P.Val < 0.05 & abs(res.lm$logFC) >= 1, 
                     ifelse(res.lm$logFC> 1,'Up', 'Down'),
                     'Stable')

write.csv(res.lm, paste0(Path, '/Data/res_limma_rpm.csv'))
