library(corrplot)


dat <- read.csv('../Ctrl.ori.csv', row.names = 1, header = T)

res <- cor(dat, method = 'pearson')

pdf('s5.corr.pdf')
corrplot(res, method = "color", tl.col ="black", order = "hclust", hclust.method = "complete", col.lim = c(0.8,1), col = colorRampPalette(c("#F8F8FF", "#F8F8FF","#FF4500"))(100), is.corr = FALSE, tl.srt = 0, addrect = 6)
dev.off()
