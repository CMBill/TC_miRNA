library(rjson)
library(tidyr)
library(miRBaseVersions.db)
library(dplyr)
Path <- getwd()

dir.create("./Data/all")
setwd(paste0(Path, '/Data/Download'))

#读取工作目录(Download)下的所有子文件夹
dir_names <- list.dirs(recursive = FALSE)
i <- 0
#读取一级子目录下文件路径，并且拷贝到父目录中all文件夹下
for (n in dir_names) {
  i <- i + 1
  x.path <- paste(n, list.files(n), sep='/')
  file.copy(x.path, '../all', recursive = FALSE)
  print(paste(i, 'in', length(dir_names)))
}

# 读取metadata
setwd(paste0(Path, '/Data'))
meta <- fromJSON(file = 'metadata.cart.2022-01-19.json')
num <- length(meta)
setwd(paste0(Path, '/Data/all'))

# 读取所有表达数据到一个数据框
exp_count <- NULL
exp_rpm <- NULL
for (i in 1:num) {
  
  # 获取一个样本，并从中得到我们需要的列
  FileName <- meta[[i]]$file_name
  EntityId <- substr(meta[[i]]$associated_entities[[1]]$entity_submitter_id, 1, 25)
  exp_of_n <- read.delim(FileName, header = TRUE)
  exp_of_n <- exp_of_n[grep('mature', exp_of_n$miRNA_region), ]
  exp_of_n <- separate(exp_of_n, into = c('State', 'ACCESSION'), col = 'miRNA_region', sep = ',')
  exp_of_n <- exp_of_n[c('miRNA_ID', 'ACCESSION', 'read_count', 'reads_per_million_miRNA_mapped')]
  colnames(exp_of_n) <- c('miRNA_ID', 'ACCESSION', 'count', 'RPM')
  
  # 从ACCESSION获取NAME
  acc <- exp_of_n$ACCESSION
  acc <- acc[!duplicated(acc)]
  ano <- miRBaseVersions.db::select(miRBaseVersions.db, keys = acc, keytype = 'VW-MIMAT-22.0', columns = c('ACCESSION', 'NAME'))
  exp_of_n <- merge(ano, exp_of_n, by = 'ACCESSION')
  
  # 计算同一NAME的miRNA表达值
  count_n <- aggregate(count ~ NAME, data = exp_of_n, sum)
  colnames(count_n) <- c('NAME', EntityId)
  rpm_n <- aggregate(RPM ~ NAME, data = exp_of_n, sum)
  colnames(rpm_n) <- c('NAME', EntityId)
  
  # 合并count
  if (is.null(exp_count)) {
    exp_count <- count_n
  }
  else{
    exp_count <- full_join(exp_count, count_n, by='NAME')
  }
  # 合并rpm
  if (is.null(exp_rpm)) {
    exp_rpm <- rpm_n
  }
  else{
    exp_rpm <- full_join(exp_rpm, rpm_n, by='NAME')
  }
}

# 将NA转化为0
exp_count[is.na(exp_count)] <- 0
exp_rpm[is.na(exp_rpm)] <- 0

write.csv(exp_count, paste0(Path, '/Data/exp_count.csv'), row.names = FALSE)
write.csv(exp_rpm, paste0(Path, '/Data/exp_rpm.csv'), row.names = FALSE)

# 以entities_id的第14-15位与10大小判断group
EntitiesId <- colnames(exp_count)
SampleGroup <- data.frame('EntitiesId'=EntitiesId)

SampleGroup$Group <- 'normal'
SampleGroup[substr(SampleGroup$EntitiesId, 14, 15) < 10, 'Group'] <- 'cancer'
SampleGroup$Group <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))
SampleGroup <- SampleGroup[-1, ]

write.csv(SampleGroup, paste0(Path, '/Data/SampleGroup.csv'), row.names = FALSE)

# 去除低表达基因，表达数小于等于0的样本超过四分之一的样本的基因被去除
# 即表达数等于0的样本小于等于四分之一的保留
selected_count <- exp_count[rowSums(exp_count <= 0) <= ncol(exp_count)/4, ]
selected_rpm <- exp_rpm[rowSums(exp_rpm <= 0) <= ncol(exp_rpm)/4, ]

write.csv(selected_count, paste0(Path, '/Data/selected_count.csv'), row.names = FALSE)
write.csv(selected_rpm, paste0(Path, '/Data/selected_rpm.csv'), row.names = FALSE)
