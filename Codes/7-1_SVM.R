library(pROC)
library(e1071)

Path <- getwd()

s.count <- read.csv('./Data/s_count.csv', check.names = FALSE, row.names = 1)
s.rpm <- read.csv('./Data/s_rpm.csv', check.names = FALSE, row.names = 1)
degs <- read.table('./Data/degs.txt')
sub <- read.table('./Data/degs.sub.txt')

SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'))
SampleGroup[, 'Group'] <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

# s.group <- sub.group[sub.group$Group != 'normal', ]
rownames(SampleGroup) <- SampleGroup$EntitiesId

# s.group$Group <- factor(s.group$Group, levels = c('PTC', 'FV-PTC'), labels = c('PTC', 'FVPTC'))

s.group <- SampleGroup[, -1 ,drop = FALSE]

# 提取两种亚型的表达样本
s.rpm <- s.rpm[, colnames(s.rpm) %in% rownames(s.group)]

rpm <- t(s.rpm)
rpmL <- rownames(rpm)
s.groupL <- rownames(s.group)

samplesL <- intersect(rpmL, s.groupL)

rpm <- rpm[samplesL, , drop = FALSE]
s.group <- s.group[samplesL, , drop = FALSE]

## 处理
# 分tran和test

# 选出degs的表达
rpm.d <- rpm[, degs$V1]

rpm1 <- cbind(rpm.d, s.group)
index <- sample(2, nrow(rpm1), replace = TRUE, prob = c(0.7, 0.3))
tran <- rpm1[index == 1,]
test <- rpm1[index == 2,]

svm.tran <- svm(Group~., data = tran, type = 'C', kernel = 'radial')
pre.svm <- predict(svm.tran, newdata = test)
obs.p.svm <- data.frame(prob = pre.svm, obs = test$Group)
con.mat.svm <- table(test$Group, pre.svm, dnn = c('R', 'P'))
svm.roc <- roc(test$Group, as.numeric(pre.svm))
plot(svm.roc, print.auc = TRUE, auc.polygon = TRUE, grid = c(0.1, 0.2), 
     grid.col = c("green", "red"), max.auc.polygon = TRUE, auc.polygon.col = "skyblue", 
     print.thres = TRUE, main = 'SVM模型ROC曲线 kernel = radial')
