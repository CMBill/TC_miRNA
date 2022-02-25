library(rjson)

Path <- getwd()
setwd(paste0(Path, '/Data'))
clinical <- fromJSON(file = 'clinical.cart.2022-01-19.json')
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
