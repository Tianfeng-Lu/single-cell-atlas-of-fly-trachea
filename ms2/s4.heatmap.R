library(pheatmap)


dat <- read.csv('all.final.pearson.csv', row.names = 1, header = T)
dat <- dat[c('PC', 'DB', 'SB', 'VB', 'TC', 'ASP', 'DT', 'LT', 'GB'),]

pdf('all.simi.pdf')
pheatmap(dat, scale = 'column', cluster_rows = F, cluster_cols = F, cellwidth = 30, cellehight = 20)
dev.off()
