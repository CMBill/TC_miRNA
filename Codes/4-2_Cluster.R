library(ConsensusClusterPlus)

Path <- getwd()

##
##

s.rpm <- read.csv(paste0(Path, '/Data/s_rpm.csv'), row.names = 1, check.names = FALSE)
SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'))
cancer.case <- SampleGroup[SampleGroup$Group == 'cancer', 'EntitiesId']
deg.rpm <- s.rpm[degs, colnames(s.rpm) %in% cancer.case]
d.rpm <- as.matrix(deg.rpm)
write.csv(deg.rpm, paste0(Path, '/Data/deg_rpm.csv'))

title <- paste0(Path, '/Data/cluster')
results <- ConsensusClusterPlus(d.rpm, maxK = 10,
                                reps = 100, pItem = 0.8,
                                pFeature = 0.8,
                                clusterAlg = "hc",
                                seed=100,
                                distance = "pearson",
                                title = title,
                                plot = "png")
icl <- calcICL(results, title = title,plot = "png")
