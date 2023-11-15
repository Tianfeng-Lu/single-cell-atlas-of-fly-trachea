## main analyses for 10x sequence data
source("./cleanRwrappers.R")

# load processed data
ctrl_filtered <- readRDS("./ctrl_filtered.rds")
HSD_filtered <- readRDS("./HSD_filtered.rds")

# UMAP for ND in Fig 1
umapplot3(ctrl_filtered)

# tracksplot for marker gene

ctrl_filtered@active.assay <- "RNA"
df <- FetchData(ctrl_filtered, c("celltype","salm","kni","ct", "vn", "wg", "bs"), slot = "data")
df$barcode <- rownames(df)
dfData <- reshape2::melt(df, id.vars = c('barcode', "celltype"), value.name="expression", variable.name='symbol')

ggplot(dfData,aes(x = barcode, y = expression, group=celltype, fill=celltype)) + ggalt::geom_horizon(colour = NA, size = 5, bandwidth = 10) + facet_grid(symbol~celltype, labeller = label_value, scales = "free", space = "free_x")  + xlab('glycosylation related proteins in HSD')  + ylab('') + theme_bw() + 
  scale_fill_manual(values = colors_list) + theme(strip.background.y = element_blank() ,strip.text.x = element_text(size = 8, color='black'), strip.text.y = element_text(hjust = 0 ,angle = 0,size = 12, color='black'), axis.text.y = element_blank(), axis.text.x = element_blank(), panel.grid = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.ticks.y = element_blank(),axis.ticks.x = element_blank(), panel.spacing.x = unit(0,"lines"), panel.spacing.y = unit(0.5,"lines"))

ctrl_filtered@active.assay <- "SCT"

# gene expression heatmap
Idents(ctrl_filtered) <- ctrl_filtered$celltype
ctrl.markers <- FindAllMarkers(ctrl_filtered, min.pct = 0.1, 
                                    logfc.threshold = 0.5, only.pos = T)
top5markers <- ctrl.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
dhm(top5markers$gene, ctrl_filtered,text_size = 6, angle=0)

# feature plot 

multi_featureplot(c("salm","kni","ct","wg","vn","bs","Tom","serp"), ctrl_filtered)

# Dot plot in Fig 3

Dotplot(c("exp","kni"), ctrl_filtered) + xlab("Gene") + ylab("Cell Type")
Dotplot(c("bs","mirr","bru2","wun"), ctrl_filtered) + xlab("Gene") + ylab("Cell Type")


# upset plot
fly_merge_filtered <- readRDS("./fly_merge_filtered.rds")

devtools::load_all("./UpSetR-master", export_all = TRUE)
Idents(ctrl_filtered) <- ctrl_filtered$celltype

all_selected_TFs <- list('ASP'=c("ken","klu","Aef1","slbo"),'DB'=c('BEAF-32','tgo',"kni","HDAC1"),'DT'=c("crol","Dref","Dp",'salm'),'GB'=c('ERR','Rel',"CG7372","dl"),'LT'=c("sima","abd-A","Usf","Hnf4"),'PC'=c('bs','pnt','ct','Trf2'),'SB'=c("wg","Atf3","grn","ct"),"TC"=c('bi',"Usf","ara","ct"),"VB"=c("Myc","Max","BEAF-32","vn"))

for(idx in seq(length(levels(Idents(ctrl_filtered))))){
  
  cell_cluster <- levels(Idents(ctrl_filtered))[idx]
  selected_TFs <- all_selected_TFs[[cell_cluster]]
  temp <- subset(ctrl_filtered, idents=cell_cluster)
  df <- FetchData(temp, selected_TFs)
  ls <- lapply(df, function(x){which(x>0)})
  no_expressed <- which(rowSums(df)==0)
  names(no_expressed) <- NULL
  ls$none <- no_expressed
  
  # svg(paste0("./figures/upset_plot/",cell_cluster,".svg"))
  print(upset(fromList(ls), order.by = "freq", keep.order = F,text.scale=2, matrix.color = "black",main.bar.color='black',point.size=4, sets.bar.color = "black",mainbar.y.label='#cell',sets.x.label='marginal size', percentage = TRUE, mainbar.y.max=NULL, nintersects = 12))
  # dev.off()
}

# split violin plot to compare HSD-ND gene expression
viotheme <- theme(plot.title = element_text(size = 17,color="black",hjust = 0.5),
                  axis.title = element_text(size = 17,color ="black"), 
                  axis.text = element_text(size = 17,color = "black"),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  axis.text.x = element_text(size= 17, angle = 0),
                  panel.grid=element_blank(),
                  legend.position = "top",
                  legend.text = element_text(size = 17),
                  legend.title= element_text(size = 17))

fly_merge_filtered@active.assay <- "SCT"
merge_expr <- FetchData(subset(fly_merge_filtered, celltype=='PC'), vars = c("fng", "kud", "mmy", "Pgant35A", "condition"), slot = "data")
merge_expr_melt <- reshape2::melt(merge_expr)
Data_summary <- Rmisc::summarySE(merge_expr_melt, measurevar="value", groupvars=c("condition","variable"))

ggplot(merge_expr_melt,aes(x = variable, y = value, fill = condition)) +
  geom_split_violin(trim=T, color="white", scale = "area") + 
  geom_point(data = Data_summary,aes(x = variable, y = value), pch=19,
             position=position_dodge(0.2), size = 1) + #绘制均值位置
  geom_errorbar(data = Data_summary, aes(ymin = value-ci, ymax = value+ci), 
                width = 0.05, 
                position= position_dodge(0.2), #误差线位置，和均值位置相匹配
                color="black",
                alpha = 0.7,
                size= 0.5) + 
  stat_compare_means(aes(group = condition),
                     label = "p.format",
                     method = "wilcox.test", size = 6,
                     label.y = max(merge_expr_melt$value),
                     hide.ns = F) + 
  scale_fill_manual(values = c("#b1d6fb", "#fd9999")) +
  labs(y=("gene expression"), x=NULL, title = "Split violin") + 
  theme_classic() + viotheme