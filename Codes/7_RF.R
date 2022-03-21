library(randomForest)

Path <- getwd()

sub.count <- read.csv('./Data/s_count.csv', check.names = FALSE, row.names = 1)
sub.rpm <- read.csv('./Data/s_rpm.csv', check.names = FALSE, row.names = 1)
sub.group <- read.csv('./Data/SampleGroup.csv')
degs <- read.table('./Data/degs.txt')
sub <- read.table('./Data/degs.sub.txt')

SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'))
SampleGroup[, 'Group'] <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

# s.group <- sub.group[sub.group$Group != 'normal', ]
rownames(SampleGroup) <- SampleGroup$EntitiesId

# s.group$Group <- factor(s.group$Group, levels = c('PTC', 'FV-PTC'), labels = c('PTC', 'FVPTC'))

s.group <- SampleGroup[, -1 ,drop = FALSE]

# 提取两种亚型的表达样本
s.rpm <- sub.rpm[, colnames(sub.rpm) %in% rownames(s.group)]

rpm <- t(s.rpm)
rpmL <- rownames(rpm)
s.groupL <- rownames(s.group)

samplesL <- intersect(rpmL, s.groupL)

rpm <- rpm[samplesL, , drop = FALSE]
s.group <- s.group[samplesL, , drop = FALSE]

# rf <- randomForest(rpm, s.group$Group)
## 患癌 
rpm.d <- rpm[, degs$V1]

rpm1 <- cbind(rpm.d, s.group)
index <- sample(2, nrow(rpm1), replace = TRUE, prob = c(0.7, 0.3))
tran <- rpm1[index == 1,]
test <- rpm1[index == 2,]

set.seed(997)

rf1 <- randomForest(tran[,1:97], tran[,98], ntree = 200,importance = T, proximity = T)
pred <- predict(rf1, newdata = test)

roc <- roc(as.ordered(test$Group), as.ordered(pred))

plot(roc, print.auc=T, auc.polygon=T, 
     max.auc.polygon=T, auc.polygon.col="white",print.thres=T)

importa <- data.frame(importance(rf1), check.names = F)

varImpPlot(rf1, n.var = min(20, nrow(rf1$importance)),
           main = 'Top 20 - variable importance')

MDSplot(rf1, tran$Group, pch = as.numeric(tran$Group))

# Fold
res.cv <- rfcv(tran[-ncol(tran)], tran$Group, cv.fold = 10)

tran.cv <- replicate(5, rfcv(tran[-ncol(tran)], tran$Group, cv.fold = 5), simplify = F)
error.cv=sapply(tran.cv,"[[","error.cv")
matplot(tran.cv[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type = "l", 
        lwd = c(2,rep(1,ncol(error.cv))), col = 1, lty = 1, log="x", 
        xlab = "Number of variables", ylab="CV Error")

## PTC FVPTC
