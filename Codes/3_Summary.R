# 不使用limma

library(tinyarray)
library(patchwork)
library(ggplot2)
library(dplyr)
Path <- getwd()

# 读取
s.count <- read.csv(paste0(Path, '/Data/s_count.csv'), row.names = 1, check.names = FALSE)
SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'), row.names = 1)
group <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

res.DESeq2 <- read.csv(paste0(Path, '/Data/res_DESeq2_rpm.csv'), row.names = 1)
res.edgeR <- read.csv(paste0(Path, '/Data/res_edgeR_rpm.csv'), row.names = 1)
res.lm <- read.csv(paste0(Path, '/Data/res_limma_rpm.csv'), row.names = 1)

res.sum <- data.frame(row.names = c('Down', 'Stable', 'Up'),
                      DESeq2 = as.integer(table(res.DESeq2$tag)),
                      edgeR = as.integer(table(res.edgeR$tag)),
                      limma = as.integer(table(res.lm$tag))
                      )

count.log <- log(s.count + 1)
count.scaled <- as.data.frame(t(scale(t(s.count))))

# 获取各自的DEGs
degs.DESeq2 <- row.names(res.DESeq2[res.DESeq2$tag != 'Stable', ])
degs.DESeq2.up <- row.names(res.DESeq2[res.DESeq2$tag == 'Up', ])
degs.DESeq2.down <- row.names(res.DESeq2[res.DESeq2$tag == 'Down', ])

degs.edgeR <- row.names(res.edgeR[res.edgeR$tag != 'Stable', ])
degs.edgeR.up <- row.names(res.edgeR[res.edgeR$tag == 'Up', ])
degs.edgeR.down <- row.names(res.edgeR[res.edgeR$tag == 'Down', ])

degs <- union(degs.DESeq2, degs.edgeR)
degs.up <- union(degs.DESeq2.up, degs.edgeR.up)
degs.down <- union(degs.DESeq2.down, degs.edgeR.down)

# 火山图
v.DESeq2 <- ggplot(res.DESeq2, aes(x = log2FoldChange, y = -log10(padj), colour = tag)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c('#4154e8', 'grey', '#e87d00')) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = 'black', lwd = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 5, col = 'black', lwd = 0.4) +
  labs(x = 'log2FC', y = 'log10FDR', title = 'DESeq2') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(1, 0.5), legend.justification = c(1, 1),
        legend.background = element_blank(), legend.title = element_blank())

v.edgeR <- ggplot(res.edgeR, aes(x = logFC, y = -log10(FDR), colour = tag)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c('#4154e8', 'grey', '#e87d00')) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = 'black', lwd = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 5, col = 'black', lwd = 0.4) +
  labs(x = 'log2FC', y = 'log10FDR', title = 'edgeR') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(1, 0.5), legend.justification = c(1, 1),
        legend.background = element_blank(), legend.title = element_blank())

# pca
pca.DESeq2 <- draw_pca(count.log[degs.DESeq2, ], group)
pca.edgeR <- draw_pca(count.log[degs.edgeR, ], group)
pca.lm <- draw_pca(count.log[degs.lm, ], group)
pca.a <- draw_pca(count.log[degs, ], group)

# 热图
h.DESeq2 <- draw_heatmap(count.log[degs.DESeq2, ], group, legend = TRUE, n_cutoff = 2)
h.edgeR <- draw_heatmap(count.log[degs.edgeR, ], group, legend = TRUE, n_cutoff = 2)
h.lm <- draw_heatmap(count.log[degs.lm, ], group, legend = TRUE, n_cutoff = 2)
h.a <- draw_heatmap(count.log[degs, ], group, legend = TRUE, n_cutoff = 2, annotation_legend = TRUE, main = 'Heat map')

# 韦恩图
up <- list(DESeq2 = degs.DESeq2.up, edgeR = degs.edgeR.up)
venn.up <- draw_venn(up, 'Up genes')
down <- list(DESeq2 = degs.DESeq2.down, edgeR = degs.edgeR.down)
venn.down <- draw_venn(down, 'Down genes')

# 拼图
v.all <- v.DESeq2 + v.edgeR
degs.all.p <- (venn.up + venn.down) / (h.a + pca.a) 

(v.DESeq2 + v.edgeR) / (venn.up + venn.down + pca.a) / h.a
