# subclustering analysis for PC DB VB TC

source("./cleanRwrappers.R")


ctrl_filtered <- readRDS("./ctrl_filtered.rds")
HSD_filtered <- readRDS("./HSD_filtered.rds")

# PC subclusters
Idents(ctrl_filtered) <- ctrl_filtered$celltype
ctrl_filtered <- FindNeighbors(ctrl_filtered, dims = 1:20)
ctrl_filtered <- FindSubCluster(ctrl_filtered, "PC", "SCT_snn", resolution = 0.15)

Idents(ctrl_filtered) <- ctrl_filtered$celltype
PC <- subset(ctrl_filtered, idents = "PC")
PC <- subset(PC, cells = CellSelector(UMAPPlot(PC))) # crop PC cells manually
Idents(PC) <- factor(PC$sub.cluster)
umapplot3(PC) + coord_fixed(ratio=1.0)+ labs(fill = "sub.cluster")

PC_sub_markers <- FindAllMarkers(PC, logfc.threshold = 0.5, min.diff.pct = 0.0, only.pos = TRUE) # pseudocount.use=0
top10markers <- PC_sub_markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)

dhm(top10markers$gene, PC)


## ORA analysis for PC subclusters
library(org.Dm.eg.db)
PC_sub_markers_GO <- FindAllMarkers(PC, logfc.threshold = log2(1.5), min.diff.pct = 0.1, pseudocount.use=0, only.pos = T) %>% filter(p_val_adj<0.01) # extended geneset for ORA 

# perform GO enrichment for each cluster
for(clus in unique(PC_sub_markers_GO$cluster)){
  genes <- bitr(filter(PC_sub_markers_GO, cluster==clus)$gene,
                fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Dm.eg.db
  )
  enrich.go <- enrichGO(
    gene = genes$ENTREZID,
    OrgDb = org.Dm.eg.db,
    keyType = "ENTREZID",
    ont = "BP", 
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    minGSSize = 10,
    readable = TRUE
  )
  
  df1 <- compute_enrichment_factor(enrich.go) %>% filter((p.adjust<0.01) & (Count>3))%>% arrange(desc(compute_ratio(enrichment_factor)))
  df2 <- df1 %>% filter(str_length(Description)<50) %>% distinct(geneID, .keep_all = TRUE) 
  
  print(ggplot(head(df2, 15), aes(x=enrichment_factor, y=factor(Description, levels = rev(Description)), fill=p.adjust)) + geom_col() + scale_fill_gradient(low = "#fd9999", high = "#b1d6fb") + 
          xlim(0,NA) + labs(x="enrichment factor", y="PC subcluster GO terms") + theme_bw()  + 
          theme(text = element_text(colour = "black", size = 16), 
                plot.title = element_text(size = 16,color="black",hjust = 0.5),
                axis.title = element_text(size = 16,color ="black"), 
                axis.text = element_text(size= 16,color = "black"), panel.grid=element_blank()) + ggtitle(clus))
}

# DT

Idents(ctrl_filtered) <- ctrl_filtered$celltype
ctrl_filtered <- FindNeighbors(ctrl_filtered, dims = 1:20)
ctrl_filtered <- FindSubCluster(ctrl_filtered, "DT", "SCT_snn", resolution = 0.15)
Idents(ctrl_filtered) <- ctrl_filtered$celltype
DT <- subset(ctrl_filtered, idents = "DT")
DT <- subset(DT, cells = CellSelector(UMAPPlot(DT))) # crop DT cells manually, remove outline cell
Idents(DT) <- factor(DT$sub.cluster)
umapplot3(DT) + coord_fixed(ratio=1.0)+ labs(fill = "sub.cluster")

# DB
Idents(ctrl_filtered) <- ctrl_filtered$celltype
ctrl_filtered <- FindNeighbors(ctrl_filtered, dims = 1:20)
ctrl_filtered <- FindSubCluster(ctrl_filtered, "DB", "SCT_snn", resolution = 0.15)
Idents(ctrl_filtered) <- ctrl_filtered$celltype
DB <- subset(ctrl_filtered, idents = "DB")
DB <- subset(DB, cells = CellSelector(UMAPPlot(DB)))
Idents(DB) <- factor(DB$sub.cluster)
umapplot3(DB) + coord_fixed(ratio=1.0)+ labs(fill = "sub.cluster")

Dotplot(c("fzr", "stg", "mei-41"), DT, scale=F)

# SB

Idents(ctrl_filtered) <- ctrl_filtered$celltype
ctrl_filtered <- FindNeighbors(ctrl_filtered, dims = 1:20)
ctrl_filtered <- FindSubCluster(ctrl_filtered, "SB", "SCT_snn", resolution = 0.15)
Idents(ctrl_filtered) <- ctrl_filtered$sub.cluster 
Idents(ctrl_filtered) <- ctrl_filtered$celltype
SB <- subset(ctrl_filtered, idents = "SB")
SB <- subset(SB, cells = CellSelector(UMAPPlot(SB))) 
Idents(SB) <- factor(SB$sub.cluster)
umapplot3(SB) + coord_fixed(ratio=1.0)+ labs(fill = "sub.cluster")

SB_sub_markers <- FindAllMarkers(SB, logfc.threshold = log2(1.2), min.diff.pct = 0.1, only.pos = FALSE)
top10markers <- SB_sub_markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)

dhm(top10markers$gene, SB, angle=0)

# VB
Idents(ctrl_filtered) <- ctrl_filtered$celltype
ctrl_filtered <- FindNeighbors(ctrl_filtered, dims = 1:20)
ctrl_filtered <- FindSubCluster(ctrl_filtered, "VB", "SCT_snn", resolution = 0.15)
Idents(ctrl_filtered) <- ctrl_filtered$sub.cluster 
Idents(ctrl_filtered) <- ctrl_filtered$celltype
VB <- subset(ctrl_filtered, idents = "VB")
VB <- subset(VB, cells = CellSelector(UMAPPlot(VB)))

Idents(VB) <- factor(VB$sub.cluster)
umapplot3(VB) + coord_fixed(ratio=1.0)+ labs(fill = "sub.cluster")
dhm(top10markers$gene, VB, angle=0)