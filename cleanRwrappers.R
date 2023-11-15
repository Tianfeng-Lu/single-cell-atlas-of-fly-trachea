# import packages
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
library(SeuratObject)
library(dplyr)

# color setting
colors_list <- c(
    "#6dc0a6", "#e2b398", "#e2a2ca", "#d1eba8", "#b1d6fb",
    "#fd9999", "#fbd69d", "#e28372", "#aeb3fb", "#c19efa", RColorBrewer::brewer.pal(9, name = "Set1")
)
aero_colors_list <- as.character(lapply(colors_list, paste0, "80"))

message(
"--------- Color Schemes ---------
For ORA  low = #fd9999, high = #b1d6fb,
For DEA  low = #1E90FF, high = #ff2121,
For cluster-specific colors: consider to construct a color mapping:
color_mapping <- colors_list[1:9]
names(color_mapping) <- ctrl_filtered$celltype %>% levels()
"
)

# visualization
umapplot3 <- function(obj, group.by = NULL){
  if (is.null(group.by)){
    group.by = "ident" # select active ident
  }
  df2 <- FetchData(obj, c("UMAP_1", "UMAP_2", group.by))
  ggobj <- ggplot(df2) + geom_point(aes(x=`UMAP_1`, y=`UMAP_2`, fill=eval(parse(text = group.by))), alpha=1, shape=21, color="white", size=2, stroke=0.2) + scale_fill_manual(values = colors_list, name = "cell type") + coord_fixed(ratio=0.8) + theme_void() + guides(fill = guide_legend(override.aes = list(size = 5)))
  LabelClusters(plot = ggobj, id = group.by, repel = F, size = 3)
  
}

featureplot <- function(gene, obj, cols = c("lightgrey", "#ff2121"), label = T, pt.size = 0.8, ...) {
    FeaturePlot(obj,
        pt.size = pt.size, features = gene, cols = cols,
        label = label, ...
    ) + theme(
        axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank()
    )
}

doublegene_featureplot <- function(gene, seuobj, ...) {
    FeaturePlot(seuobj,
        pt.size = 0.8, features = gene,
        cols = c("lightgrey", "deeppink", "#1E90FF"), label = TRUE, blend = T, ...
    ) + theme(
        axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()
    )
}


multi_featureplot <- function(genes, seuobj,
                              nrow = ceiling(sqrt(length(genes))),
                              ncol = ceiling(sqrt(length(genes))),
                              labels = "AUTO", ...) {
    genes <- intersect(genes, rownames(seuobj))
    cowplot::plot_grid(
        plotlist = lapply(genes, featureplot, seuobj, ...),
        labels = labels, ncol = ncol, nrow = ncol, label_size = 18
    )
}

Dotplot <- function(gene, seuobj, ...) {
    DotPlot(seuobj, features = gene, cols = c("lightgrey", "#ff2121"), col.min = 0, ...) +
         theme_classic() + theme(text = element_text(colour = "black", size = 16), 
                                 plot.title = element_text(size = 16,color="black",hjust = 0.5),
                                 axis.title = element_text(size = 16,color ="black"), 
                                 axis.text = element_text(size = 16,color = "black"))
}

violin_plot <- function(gene, seuobj, ...) {
    VlnPlot(seuobj,
        features = gene, stack = T, pt.size = 0, cols = colors_list, 
        direction = "horizontal", 
        x.lab = "", y.lab = "", ...) 
        theme(
            axis.text.x = element_blank(), axis.text = element_text(size= 16,color = "black"),
            axis.ticks.x = element_blank()
        )+ 
       theme(text = element_text(colour = "black", size = 18))
}

dhm <- function(gene, seuobj,text_size = 6, ...) {
    DoHeatmap(seuobj, features = gene, group.colors = colors_list, size = text_size, ...) +
        scale_fill_gradientn(colors = c("#1E90FF", "white", "#ff2121")) + theme(axis.text = element_text(colour = "black", size = 12))
}


dhm2 <- function(gene, seuobj, genes_to_show, session = "session", ...){
  
  mat <- GetAssayData(seuobj, slot = "scale.data")
  
  cluster_info <- sort(Idents(seuobj))
  heatmapgenes <- intersect(gene,rownames(mat))
  mat <- as.matrix(mat[heatmapgenes,names(cluster_info)])
  gene_pos <- match(genes_to_show, rownames(mat))
  row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,labels = genes_to_show))
  
  col <- colors_list
  names(col) <- levels(cluster_info)
  top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=col),labels = levels(cluster_info),
                                                   labels_gp = gpar(cex=1,col='black'))) 

  col_fun <-  colorRamp2(c(-2, 1, 4), c("#1E90FF", "white", "#ff2121"))
  
  Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, 
          show_column_names = FALSE, show_row_names = FALSE,
          column_split = cluster_info, top_annotation = top_anno,  
          column_title = NULL, right_annotation = row_anno, 
          heatmap_legend_param = list(
            title='Expression', title_position='leftcenter-rot'), col = col_fun)
}

# seurat4 function
FindSubCluster <- function(object,
                           cluster,
                           graph.name = "RNA_snn",
                           subcluster.name = "sub.cluster",
                           resolution = 0.3,
                           algorithm = 1) {
    sub.cell <- WhichCells(object = object, idents = cluster)
    sub.graph <- as.Graph(x = object[[graph.name]][sub.cell, sub.cell])
    sub.clusters <- FindClusters(
        object = sub.graph,
        resolution = resolution,
        algorithm = algorithm
    )
    sub.clusters[, 1] <- paste(cluster, sub.clusters[, 1], sep = "_")
    object[[subcluster.name]] <- as.character(x = Idents(object = object))
    object[[subcluster.name]][sub.cell, ] <- sub.clusters[, 1]
    return(object)
}

# average expression matrix each cluster
get_averexpr_mat_cluster <- function(seuobj, type = "data") {
    aver_expr_mat <- c()
    for (cluster_name in levels(seuobj@active.ident))
    {
        temp <- as.matrix(GetAssayData(seuobj, slot = type)[, WhichCells(seuobj, idents = cluster_name)])
        aver_expr_mat <- cbind(aver_expr_mat, apply(temp, 1, mean))
    }
    colnames(aver_expr_mat) <- levels(seuobj@active.ident)
    return(aver_expr_mat)
}

# return highvar genes countmatrix/datamatrix/cpmmatrix
get_data_table <- function(seuobj, highvar = T, type = "counts") {
    if (type == "cpm") {
        countsMat <- as.matrix(GetAssayData(seuobj, slot = "counts"))
    } else {
        countsMat <- as.matrix(GetAssayData(seuobj, slot = type))
    }
    if (highvar == T) {
        # genelist <- GetAssayData(seuobj, slot = "var.features") 
        genelist <- seuobj@assays[["SCT"]]@var.features
        countsMat <- countsMat[genelist, ] 
    }
    cpm <- apply(countsMat, 2, function(x) {
        x / sum(x) * 1000000
    })
    if (type == "cpm") {
        return(cpm)
    } else {
        return(countsMat)
    }
}


scatterplot <- function(gene, subseuobj, type = "sub.cluster") {
    embedding <- FetchData(object = subseuobj, vars = c("UMAP_1", "UMAP_2", type, gene))
    embedding$barcode <- paste(rownames(embedding))
    plot_ly(x = embedding$UMAP_1, y = embedding$UMAP_2, z = embedding$gene, type = "scatter3d", size = 1, color = embedding$sub.cluster, colors = colors_list, mode = "markers")
}



sankey_plot <- function(confuse_matrix, label1 = dimnames(confuse_matrix)$pre, label2 = dimnames(confuse_matrix)$true, session = "session") {
    sources <- rep(0:(length(label1) - 1), each = length(label2)) # 注意这里的each和times的区别
    colors <- rep(aero_colors_list[1:length(label1)], each = length(label2))
    targets <- rep(length(label1) + 0:(length(label2) - 1), times = length(label1))

    plot_ly(
        type = "sankey", orientation = "h",
        node = list(
            # label = c(label1, label2),
          label = NULL,
            color = c(colors_list[1:length(label1)],colors_list[1:length(label2)]), pad = 15, thickness = 30,
            line = list(color = "black", width = 1)
        ),
        link = list(
            source = sources, target = targets,
            value = as.numeric(confuse_matrix),
            color = colors
        )
    ) %>% layout(title = session, font = list(family = "Arial", size = 20, color = "black"))
}


confuse_bubblemat <- function(confuse_matrix_prop, label1, label2, session = "session") {
    prop <- as.numeric(confuse_matrix_prop)
    data <- expand.grid(x = label1, y = label2) %>% bind_cols(prop = prop)
    plot <- ggplot(data, aes(x = x, y = y, colour = prop, size = prop)) +
        geom_point() +
        scale_size_continuous(range = c(0, 10)) +
        labs(x = "clusters", y = "inferred from") +
        theme_bw()

    ggsave(paste0(session, ".svg"), plot = plot, device = svg, width = 5, height = 4)
}



mytheme <- theme(plot.title = element_text(size = 15,color="black",hjust = 0.5),
                 axis.title = element_text(size = 15,color ="black"), 
                 axis.text = element_text(size = 15,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 15),
                 legend.title= element_text(size= 15)) 

mytheme2 <- theme(plot.title = element_text(size = 15,color="black",hjust = 0.5),
                 axis.title = element_text(size = 15,color ="black"), 
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(size = 15,color = "black"),
                 panel.background = element_rect(fill = "white"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 panel.grid=element_blank(),
                 legend.position = "bottom",
                 legend.text = element_text(size= 15),
                 legend.title= element_text(size= 15)) 

# https://stackoverflow.com/a/45614547
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

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

scatter_theme <- theme(axis.title = element_text(size = 20,color = "black"),
                       axis.text = element_text(size = 20,color = "black"),
                       axis.line = element_line(size = 1),
                       axis.ticks = element_line(size = 1),
                       legend.key.size = unit(2,"cm"),
                       title = element_text(size = 20))

ridgetheme <- theme(plot.title = element_text(size = 15,color="black",hjust = 0.5),
                    axis.title = element_text(size = 15,color ="black"), 
                    axis.text = element_text(size = 15,color = "black"),
                    panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    axis.text.x = element_text(angle = 0, hjust = 1),
                    panel.grid=element_blank(),
                    legend.position = "top",
                    legend.text = element_text(size= 15),
                    legend.title= element_text(size= 15))

bartheme <- theme(plot.title = element_text(size = 18, color="black",hjust = 0.5),
                  axis.title = element_text(size = 18,color ="black"),
                  axis.text = element_text(size = 18,color = "black"),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1 ),
                  panel.grid=element_blank(),
                  legend.position = "top",
                  legend.text = element_text(size = 18),
                  legend.title= element_text(size = 18))


reg_vector <- function(x){
  paste0("c(",toString(x),")")
}


# for ORA -----------------------------------------------------------------

compute_ratio <- function(x){sapply(x, function(y){parse(text=y) %>% eval()})}
compute_enrichment_factor <- function(x) {return(x@result %>% mutate(enrichment_factor=compute_ratio(GeneRatio)/compute_ratio(BgRatio)) %>%  arrange(desc(enrichment_factor)))}