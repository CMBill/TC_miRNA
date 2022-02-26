library(rjson)

Path <- getwd()
setwd(paste0(Path, '/Data'))

clinical <- fromJSON(file = 'clinical.cart.2022-01-19.json')
SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'))
s.rpm <- read.csv(paste0(Path, '/Data/s_rpm.csv'), row.names = 1, check.names = FALSE)
s.count <- read.csv(paste0(Path, '/Data/s_count.csv'), row.names = 1, check.names = FALSE)

num <- length(clinical)
clin <- data.frame('case' = character(), 'subtype' = character())

# 从临床数据读取亚型
for (i in 1:num) {
  type <- clinical[[i]]$diagnoses[[1]]$primary_diagnosis
  case <- substr(clinical[[i]]$diagnoses[[1]]$submitter_id, 1, 12)
  clin[i, ]$case <- case
  clin[i, ]$subtype <- type
}

# 提取出需要的两种亚型
s1 <- clin[clin$subtype == 'Papillary adenocarcinoma, NOS', ]
s2 <- clin[clin$subtype == 'Papillary carcinoma, follicular variant', ]
subtype <- rbind(s1, s2)

SampleGroup$case <- substr(SampleGroup$EntitiesId, 1, 12)
s.groups <- merge(SampleGroup, subtype, by = 'case')
s.groups[s.groups$Group == 'normal', 'subtype'] <- 'normal'
s.groups[s.groups$subtype == 'Papillary adenocarcinoma, NOS', 'Group'] <- 'PTC'
s.groups[s.groups$subtype == 'Papillary carcinoma, follicular variant', 'Group'] <- 'FV-PTC'
SubGroup <- s.groups[c('EntitiesId', 'Group')]

sub.rpm <- s.rpm[, colnames(s.rpm) %in% SubGroup$EntitiesId]
sub.count <- s.count[, colnames(s.count) %in% SubGroup$EntitiesId]

write.csv(SubGroup, paste0(Path, '/Data/sub_group.csv'), row.names = FALSE)
write.csv(sub.rpm, paste0(Path, '/Data/sub_rpm.csv'), row.names = TRUE)
write.csv(sub.count, paste0(Path, '/Data/sub_count.csv'), row.names = TRUE)
