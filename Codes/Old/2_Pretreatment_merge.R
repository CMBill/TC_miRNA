library(rjson)

Path <- getwd() # [1] "E:/wkh/Codes/Projects/TC_Demo"
# 读取metadata
setwd(paste0(Path, '/Data'))
meta <- fromJSON(file = 'metadata.cart.2022-01-11.json')
num <- length(meta)
setwd(paste0(Path, '/Data/all'))
exp_merge <- NULL

# 读取所有表达数据到一个数据框
for (i in 1:num) {
  FileName <- substr(meta[[i]]$file_name, 1, 49)
  EntityId <- substr(meta[[i]]$associated_entities[[1]]$entity_submitter_id, 1, 25)
  exp_of_n <- read.delim(FileName, col.names = c('ID', EntityId), header=FALSE)
  if (is.null(exp_merge)) {
    exp_merge <- exp_of_n
  }
  else{
    exp_merge <- merge(exp_merge, exp_of_n, by='ID')
  }
}

exp_data <- exp_merge[- (1:5), ]
rownames(exp_data) <- gsub("\\.(\\.?\\d*)", '', exp_data$ID)
exp_data <- exp_data[, -1]

write.csv(exp_data, paste0(Path, '/Data/Rexp_data.csv'))

# 以entities_id的第14-15位与10大小判断group
EntitiesId <- colnames(exp_data)
SampleGroup <- data.frame('EntitiesId'=EntitiesId)
SampleGroup[, 'Group'] <- 'normal'
SampleGroup[substr(SampleGroup$EntitiesId, 14, 15) < 10, 'Group'] <- 'cancer'
SampleGroup[, 'Group'] <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

write.csv(SampleGroup, paste0(Path, '/Data/RSampleGroup.csv'))

# 去除低表达基因，表达数小于10的样本超过四分之一的样本的基因被去除
# 即表达数小于10的样本小于等于四分之一的保留
exp_select <- exp_data[rowSums(exp_data<10)<=ncol(exp_data)/4, ]

write.csv(exp_select, paste0(Path, '/Data/Rexp_select.csv'))