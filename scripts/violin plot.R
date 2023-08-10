# ctrl_filtered <- ScaleData(ctrl_filtered,vars.to.regress = 'percent.mt')
# ctrl_filtered <- FindVariableFeatures(ctrl_filtered, selection.method = "vst", nfeatures = 2000)
# ctrl_filtered <- ScaleData(ctrl_filtered, features = all.genes)

ctrldata <- get_data_table(subset(ctrl_filtered, ident="DB"),type = "data",highvar = F)
HSDdata <- get_data_table(subset(HSD_filtered, ident="DB"),type = "data",highvar = F)
genes_to_show <- c('bs', 'kni')

func1 <- function(gene, sample, datable){
  data.frame(expr = datable[gene,], sample = sample, gene = gene)
}

merge_expr <- data.frame()
for (i in lapply(genes_to_show, func1,"ctrldata",ctrldata))
{
  merge_expr <- rbind(merge_expr,i)
}
for (i in lapply(genes_to_show, func1,"HSDdata",HSDdata))
{
  merge_expr <- rbind(merge_expr,i)
}

rownames(merge_expr) <- NULL
Data_summary <- Rmisc::summarySE(merge_expr, measurevar="expr", groupvars=c("sample","gene"))
# head(Data_summary)

viotheme <- theme(plot.title = element_text(size = 17,color="black",hjust = 0.5),
                  axis.title = element_text(size = 17,color ="black"), 
                  axis.text = element_text(size = 17,color = "black"),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  axis.text.x = element_text(size= 17, angle = 0),
                  panel.grid=element_blank(),
                  legend.position = "top",
                  legend.text = element_text(size= 17),
                  legend.title= element_text(size= 17))

stat_theme <- stat_compare_means(aes(group = sample),
                                 label = "p.format",
                                 method = "wilcox.test", size = 6,
                                 label.y = max(merge_expr$expr),
                                 hide.ns = F)

ggobj <- ggplot(merge_expr,aes(x = gene, y = expr, fill = sample)) +
  geom_split_violin(trim= F, color="white", scale = "area") + 
  geom_point(data = Data_summary,aes(x = gene, y= expr), pch=19,
             position=position_dodge(0.2),size= 1) + #绘制均值位置
  geom_errorbar(data = Data_summary, aes(ymin = expr-ci, ymax= expr+ci), 
                width= 0.05, 
                position= position_dodge(0.2), #误差线位置，和均值位置相匹配
                color="black",
                alpha = 0.7,
                size= 0.5) +
  scale_fill_manual(values = c("#b1d6fb", "#fd9999"))+ 
  labs(y=("gene expression"),x=NULL, title = "DB ctrl vs HSD") + 
  theme_classic()+ viotheme + stat_theme
ggobj
ggsave("./figures/DB ctrl vs HSD1.tiff", device = tiff, plot = ggobj, height = 5, width = 7)
