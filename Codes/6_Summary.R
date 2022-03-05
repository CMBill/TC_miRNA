library(tinyarray)
library(patchwork)
library(ggplot2)
library(dplyr)
Path <- getwd()
# Neg & Pos
s.count <- read.csv(paste0(Path, '/Data/s_count.csv'), row.names = 1, check.names = FALSE)
SampleGroup <- read.csv(paste0(Path, './Data/SampleGroup.csv'), row.names = 1)
group <- factor(SampleGroup[, 'Group'], levels = c('normal', 'cancer'), labels = c('normal', 'cancer'))

res.DESeq2 <- read.csv(paste0(Path, '/Data/res_DESeq2_rpm.csv'), row.names = 1)
res.edgeR <- read.csv(paste0(Path, '/Data/res_edgeR_rpm.csv'), row.names = 1)

res.sum <- data.frame(row.names = c('Down', 'Stable', 'Up'),
                      DESeq2 = as.integer(table(res.DESeq2$tag)),
                      edgeR = as.integer(table(res.edgeR$tag))
)

count.log <- log(s.count + 1)
count.scaled <- as.data.frame(t(scale(t(s.count))))

degs.DESeq2 <- row.names(res.DESeq2[res.DESeq2$tag != 'Stable', ])
degs.DESeq2.up <- row.names(res.DESeq2[res.DESeq2$tag == 'Up', ])
degs.DESeq2.down <- row.names(res.DESeq2[res.DESeq2$tag == 'Down', ])

degs.edgeR <- row.names(res.edgeR[res.edgeR$tag != 'Stable', ])
degs.edgeR.up <- row.names(res.edgeR[res.edgeR$tag == 'Up', ])
degs.edgeR.down <- row.names(res.edgeR[res.edgeR$tag == 'Down', ])

degs <- union(degs.DESeq2, degs.edgeR)
degs.up <- union(degs.DESeq2.up, degs.edgeR.up)
degs.down <- union(degs.DESeq2.down, degs.edgeR.down)

res.DESeq2.order <- res.DESeq2[order(abs(res.DESeq2$padj), decreasing = T),]
res.edgeR.order <- res.edgeR[order(abs(res.edgeR$FDR), decreasing = T),]

degs.DESeq2.top10 <- rownames(res.DESeq2.order[1:10,])
degs.edgeR.top10 <- rownames(res.edgeR.order[1:10,])

write(degs.DESeq2, './Data/degs.DESeq2.txt')
write(degs.DESeq2.up, './Data/degs.DESeq2.up.txt')
write(degs.DESeq2.down, './Data/degs.DESeq2.down.txt')
write(degs.edgeR, './Data/degs.edgeR.txt')
write(degs.edgeR.up, './Data/degs.edgeR.up.txt')
write(degs.edgeR.down, './Data/degs.edgeRdown.txt')
write(degs, './Data/degs.txt')
write(degs.up, './Data/degs.up.txt')
write(degs.down, './Data/degs.down.txt')

# PTC & FVPTC
sub.group <- read.csv('./Data/sub_group.csv')
s_ptc <- sub.group[sub.group$Group == 'PTC', ]
s_fvptc <- sub.group[sub.group$Group == 'FV-PTC', ]
s_nor <- sub.group[sub.group$Group == 'normal', ]

s.group <- rbind(s_ptc, s_fvptc)
row.names(s.group) <- s.group$EntitiesId
s.group <- s.group[, -1, drop = FALSE]

sub.count <- s.count[, colnames(s.count) %in% row.names(s.group)]
sub.count.log <- log(sub.count + 1)
sub.count.scaled <- as.data.frame(t(scale(t(sub.count))))

s_group <- s.group[colnames(sub.count), , drop = FALSE]
sgroup <- factor(s_group$Group, levels = c('PTC', 'FV-PTC'), labels = c('PTC', 'FVPTC'))

res.DESeq2.sub <- read.csv(paste0(Path, '/Data/res_DESeq2_sub.csv'), row.names = 1)
res.edgeR.sub <- read.csv(paste0(Path, '/Data/res.edgeR_sub.csv'), row.names = 1)

degs.DESeq2.sub <- row.names(res.DESeq2.sub[res.DESeq2.sub$tag != 'Stable', ])
degs.DESeq2.up.sub <- row.names(res.DESeq2.sub[res.DESeq2.sub$tag == 'Up', ])
degs.DESeq2.down.sub <- row.names(res.DESeq2.sub[res.DESeq2.sub$tag == 'Down', ])

degs.edgeR.sub <- row.names(res.edgeR.sub[res.edgeR.sub$tag != 'Stable', ])
degs.edgeR.up.sub <- row.names(res.edgeR.sub[res.edgeR.sub$tag == 'Up', ])
degs.edgeR.down.sub <- row.names(res.edgeR.sub[res.edgeR.sub$tag == 'Down', ])

degs.sub <- union(degs.DESeq2.sub, degs.edgeR.sub)
degs.up.sub <- union(degs.DESeq2.up.sub, degs.edgeR.up.sub)
degs.down.sub <- union(degs.DESeq2.down.sub, degs.edgeR.down.sub)

res.sum.sub <- data.frame(row.names = c('Down', 'Stable', 'Up'),
                      DESeq2 = as.integer(table(res.DESeq2.sub$tag)),
                      edgeR = as.integer(table(res.edgeR.sub$tag))
)

write(degs.DESeq2.sub, './Data/degs.DESeq2.sub.txt')
write(degs.DESeq2.up.sub, './Data/degs.DESeq2.up.sub.txt')
write(degs.DESeq2.down.sub, './Data/degs.DESeq2.down.sub.txt')
write(degs.edgeR.sub, './Data/degs.edgeR.sub.txt')
write(degs.edgeR.up.sub, './Data/degs.edgeR.up.sub.txt')
write(degs.edgeR.down.sub, './Data/degs.edgeRdown.sub.txt')
write(degs.sub, './Data/degs.sub.txt')
write(degs.up.sub, './Data/degs.up.sub.txt')
write(degs.down.sub, './Data/degs.down.sub.txt')

# PCA
pca.DESeq2 <- draw_pca(count.log[degs.DESeq2, ], group)
pca.edgeR <- draw_pca(count.log[degs.edgeR, ], group)
pca.a <- draw_pca(count.log[degs, ], group)

pca.DESeq2.sub <- draw_pca(sub.count.log[degs.DESeq2.sub, ], sgroup)
pca.edgeR.sub <- draw_pca(sub.count.log[degs.edgeR.sub, ], sgroup)
pca.a.sub <- draw_pca(sub.count.log[degs.sub, ], sgroup)

pca <- draw_pca(s.count, group)
pca.sub <- draw_pca(sub.count, sgroup)

# Vol
v.DESeq2 <- ggplot(res.DESeq2, aes(x = log2FoldChange, y = -log10(padj), colour = tag)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c('#0067a7', 'grey', '#ef8026')) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = 'black', lwd = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 5, col = 'black', lwd = 0.4) +
  labs(x = 'log2FC', y = 'log10FDR', title = 'DESeq2') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(1, 0.5), legend.justification = c(1, 1),
        legend.background = element_blank(), legend.title = element_blank())

v.edgeR <- ggplot(res.edgeR, aes(x = logFC, y = -log10(FDR), colour = tag)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c('#0067a7', 'grey', '#ef8026')) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = 'black', lwd = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 5, col = 'black', lwd = 0.4) +
  labs(x = 'log2FC', y = 'log10FDR', title = 'edgeR') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(1, 0.5), legend.justification = c(1, 1),
        legend.background = element_blank(), legend.title = element_blank())

v.DESeq2.sub <- ggplot(res.DESeq2.sub, aes(x = log2FoldChange, y = -log10(padj), colour = tag)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c('#0067a7', 'grey', '#ef8026')) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = 'black', lwd = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 5, col = 'black', lwd = 0.4) +
  labs(x = 'log2FC', y = 'log10FDR', title = 'DESeq2') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(1, 0.5), legend.justification = c(1, 1),
        legend.background = element_blank(), legend.title = element_blank())

v.edgeR.sub <- ggplot(res.edgeR.sub, aes(x = logFC, y = -log10(FDR), colour = tag)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c('#0067a7', 'grey', '#ef8026')) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = 'black', lwd = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 5, col = 'black', lwd = 0.4) +
  labs(x = 'log2FC', y = 'log10FDR', title = 'edgeR') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(1, 0.5), legend.justification = c(1, 1),
        legend.background = element_blank(), legend.title = element_blank())

# Venn
up <- list(DESeq2 = degs.DESeq2.up, edgeR = degs.edgeR.up)
venn.up <- draw_venn(up, 'Up genes')
down <- list(DESeq2 = degs.DESeq2.down, edgeR = degs.edgeR.down)
venn.down <- draw_venn(down, 'Down genes')

up.sub <- list(DESeq2 = degs.DESeq2.up.sub, edgeR = degs.edgeR.up.sub)
venn.up.sub <- draw_venn(up.sub, 'Up genes')
down.sub <- list(DESeq2 = degs.DESeq2.down.sub, edgeR = degs.edgeR.down.sub)
venn.down.sub <- draw_venn(down.sub, 'Down genes')


