#整理第一次复制出来的tar.gz格式表达文件
setwd("E:/wkh/Codes/Projects/TC_Demo/Data/Download")
dir.create("../all")
#读取工作目录(Download)下的所有子文件夹
dir_names <- list.dirs(recursive = FALSE)
i <- 0
#读取一级子目录下文件路径，并且拷贝到父目录中all文件夹下
for (n in dir_names) {
  i <- i + 1
  x.path <- paste(n, list.files(n), sep='/')
  file.copy(x.path, '../all', recursive = FALSE)
  print(paste(i, 'in', length(dir_names), '\r'))
}

#在all文件夹下解压缩gz格式文件并删除原文件
setwd("E:/wkh/Codes/Projects/TC_Demo/Data/all")
library(R.utils)
file.remove("annotations.txt")
file_names <- list.files(recursive = FALSE)
j <- 0
for (n in file_names) {
  j <- j + 1
  gunzip(n, remove = 'TRUE')
  print(paste(j, 'in', length(file_names), '\r'))
}
print('Unpack Done')
