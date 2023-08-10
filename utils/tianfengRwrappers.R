###
# f -- featureplot for single gene
# multi_featureplot
# dhm -- doheatmap
# doublegene_featureplot
# Dotplot
# violin_plot
# umapplot
# make_bins make pseudotime bins
# make_pseudotime_heatmap_bycells
# make_pseudotime_heatmap_average
# ridgeline_plot
# FindSubCluster
# get_averexpr_mat_cluster
# get_data_table
# GO_dotplot
# 3D_featureplot
###



library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(MySeuratWrappers)
library(ggpubr)
library(cowplot)
library(pheatmap)
library(viridis)
library(reshape2)
library(ggridges)
library(hrbrthemes)
library(clusterProfiler)
library(extrafont)
library(plotly)
library(listviewer)
library(stringr)
library(ggprism)
# library(Mfuzz)
library(IDPmisc)
library(ggrepel)
# library(future)
# plan("multiprocess", workers = 6)
library(circlize)
library(ComplexHeatmap)
library(ggnewscale)
# library(org.Hs.eg.db)
# schema(jsonedit = interactive())


colors_list <- c(
    "#6dc0a6", "#e2b398", "#e2a2ca", "#d1eba8", "#b1d6fb",
    "#fd9999", "#fbd69d", "#e28372", "#aeb3fb", "#c19efa", brewer.pal(9, name = "Set1")
)
aero_colors_list <- as.character(lapply(colors_list, paste0, "80")) # 透明化颜色

message(
"--------- Color Schemes ---------
For ORA  low = #fd9999, high = #b1d6fb,
For DEA  low = #1E90FF, high = #ff2121,
For cluster-specific colors: consider to construct a color mapping:
color_mapping <- colors_list[1:9]
names(color_mapping) <- ctrl_filtered$celltype %>% levels()
"
)

umapplot <- function(obj, group.by = NULL, label = TRUE, pt.size = 0.8, ...) {
  DimPlot(obj,
          reduction = "umap", cols = colors_list, label = label,
          pt.size = pt.size, group.by = group.by, ...) +
    xlab("UMAP 1") + ylab("UMAP 2") + guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme(axis.line = element_line(arrow = arrow(length = unit(0.2, "cm")))) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL) + scale_color_manual(values = colors_list, drop = F) #drop 参数用来绘制所有因子水平的颜色
}


umapplot2 <- function(obj, group.by = NULL, label = TRUE, ...) {
  umapplot(obj,group.by = "celltype") + theme(axis.line = element_blank(),axis.title = element_blank(),legend.position = 'none')+coord_fixed(ratio = 1) 
  
}

f <- function(gene, obj, cols = c("lightgrey", "#ff2121"), label = T, pt.size = 0.8, ...) {
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

## update 更新了feature的选取，避免找不到对应的基因报错
## labels can be "AUTO" or NA
multi_featureplot <- function(genes, seuobj,
                              nrow = ceiling(sqrt(length(genes))),
                              ncol = ceiling(sqrt(length(genes))),
                              labels = "AUTO", ...) {
    genes <- intersect(genes, rownames(seuobj))
    plot_grid(
        plotlist = lapply(genes, f, seuobj, ...),
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
        features = gene, stack = T, pt.size = 0, cols = colors_list, # 颜色
        direction = "horizontal", # horizontal 水平作图  vertical
        x.lab = "", y.lab = "", ...) + # 横纵轴不标记任何东西
        theme(
            axis.text.x = element_blank(), axis.text = element_text(size= 16,color = "black"),
            axis.ticks.x = element_blank()
        )+ # 不显示坐标刻度  
       theme(text = element_text(colour = "black", size = 18))
}

dhm <- function(gene, seuobj,text_size = 6, ...) {
    DoHeatmap(seuobj, features = gene, group.colors = colors_list, size = text_size, ...) +
        scale_fill_gradientn(colors = c("#1E90FF", "white", "#ff2121")) + theme(axis.text = element_text(colour = "black", size = 12))
} # 修正的doheatmap绘图颜色


dhm2 <- function(gene, seuobj, genes_to_show, session = "session", ...){
  
  mat <- GetAssayData(seuobj, slot = "scale.data")
  
  ##获得基因和细胞聚类信息
  cluster_info <- sort(Idents(seuobj))
  heatmapgenes <- intersect(gene,rownames(mat))
  
  ##筛选矩阵为要画热图的基因
  mat <- as.matrix(mat[heatmapgenes,names(cluster_info)])
  
  #获得要展示的基因在热图中的位置信息
  gene_pos <- match(genes_to_show, rownames(mat))
  row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,labels = genes_to_show))
  
  col <- colors_list
  names(col) <- levels(cluster_info)
  top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=col),labels = levels(cluster_info),
                                                   labels_gp = gpar(cex=1,col='black'))) ## 顶端的cluster注释
  

  col_fun <-  colorRamp2(c(-2, 1, 4), c("#1E90FF", "white", "#ff2121")) #颜色
  
  Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, 
          show_column_names = FALSE, show_row_names = FALSE,
          column_split = cluster_info, top_annotation = top_anno,  
          column_title = NULL, right_annotation = row_anno, 
          heatmap_legend_param = list(
            title='Expression', title_position='leftcenter-rot'), col = col_fun)
}


# return seuobj with pseudotime bins
make_bins <- function(seuobj, nbins = 11) {
    ordered.cells <- BiocGenerics::rank(seuobj$pseudotime)
    bins <- c()
    len <- length(seuobj$pseudotime)
    for (i in c(1:len)) {
        bins[i] <- floor(ordered.cells[i] * nbins / len)
        # 创建一个长度为细胞数量的向量，存储每个细胞分别属于的拟时bin信息
    }
    Idents(seuobj) <- as.factor(bins)
    seuobj$bins <- bins
    return(seuobj)
}

# 按照拟时分群绘制每个细胞中的基因表达
make_pseudotime_heatmap_bycells <- function(seuobj, top_specific_marker_ids, topn_markers = 40) {
    dhm(head(top_specific_marker_ids, topn_markers), seuobj, group.by = "bins")
}

# 根据拟时平均画图
# 提取表达矩阵
make_pseudotime_heatmap_average <- function(seuobj,
                                            top_specific_marker_ids,
                                            topn_markers = 40,
                                            nbins = 11) {
    mat <- c()
    # 上述方法分出来的bins是从0开始的
    for (j in c(0:(nbins - 1)))
    {
        temp <- as.matrix(GetAssayData(seuobj, slot = "scale.data")[, WhichCells(seuobj, idents = c(j))])
        mat <- cbind(mat, apply(temp, 1, mean))
    }
    colnames(mat) <- as.factor(c(0:(nbins - 1)))

    pheatmap(mat[head(top_specific_marker_ids, topn_markers), ],
        breaks = unique(c(seq(-3, 3, length = 400))),
        color = colorRampPalette(c("#1E90FF", "white", "#ff2121"))(400),
        border_color = NA, cluster_rows = FALSE, cluster_cols = FALSE,
        main = "Pseudotime Bins", angle_col = 45, show_rownames = T
    )
}


# 绘制ridgeline plot
# 首先设置阈值来二值化对特定基因在每个细胞中是否表达
# 然后观察表达这个基因的细胞的拟时值分布

ridgeline_plot <- function(genelist, seuobj, expr_threshold = 1) {
    scale_data <- GetAssayData(seuobj, slot = "scale.data")
    df <- c()
    df$pseudotime <- NULL
    df$gene <- NULL
    for (gene in genelist) {
        candidate <- rownames(as.data.frame(scale_data[gene, ][scale_data[gene, ] > expr_threshold]))
        temp <- as.data.frame(seuobj$pseudotime[candidate])
        colnames(temp) <- "pseudotime"
        rownames(temp) <- NULL
        temp$gene <- gene
        df <- rbind(df, temp)
    }
    plot1 <- ggplot(df, aes(x = pseudotime, y = gene, fill = gene)) +
        ggridges::geom_density_ridges(alpha = 0.5, show.legend = T) +
        labs(title = "DEGs with pseudotime") +
        theme_ipsum() +
        theme(
            legend.position = "none", panel.spacing = unit(0.1, "lines"),
            strip.text.x = element_text(size = 8)
        )
    plot2 <- ggplot(df, aes(x = pseudotime, y = gene, fill = ..x..)) +
        geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, alpha = 0.5) +
        scale_fill_viridis(name = seuobj$pseudotime, option = "C") +
        labs(title = "DEGs with pseudotime") +
        theme_ipsum() +
        theme(
            legend.position = "none",
            panel.spacing = unit(0.1, "lines"),
            strip.text.x = element_text(size = 8)
        )
    return(list(plot1, plot2))
}


genewithpseudotime <- function(gene, seuobj) {
    datMat <- as.matrix(GetAssayData(seuobj, slot = "data"))
    df <- data.frame(t(seuobj$pseudotime))
    df <- do.call(data.frame, lapply(df, function(x) {
        x <- x / 20
        replace(x, is.infinite(x), 0)
    }))
    rownames(df) <- "pseudotime"
    colnames(df) <- colnames(datMat)
    datMat <- rbind(df, datMat)

    test <- CreateSeuratObject(counts = datMat)
    test@reductions[["umap"]] <- seuobj@reductions[["umap"]]
    test@active.ident <- seuobj@active.ident
    # SetAssayData(seuobj,slot = "data", new.data = datMat)

    gene_name <- c(gene, "pseudotime")
    FeaturePlot(test,
        pt.size = 0.8, slot = "counts", min.cutoff = 0,
        max.cutoff = 4, features = gene_name, cols = c("lightgrey", "#1E90FF", "#ff2121"),
        label = TRUE, blend = T, blend.threshold = 0
    ) +
        theme(
            axis.line = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), axis.title = element_blank()
        )
}


# 组合pseudotime 波浪图
multi_ridgeplot <- function(genes, seuobj, xaxis = "pseudotime",
                            nrow = ceiling(sqrt(length(genes))),
                            ncol = ceiling(sqrt(length(genes))),
                            labels = "AUTO", ...) {
  df <- FetchData(seuobj, vars = c("pseudotime",genes))
  data <- cbind(df,index = 1:nrow(df), cluster = seuobj$celltype) %>% na.omit()
  genes <- intersect(genes, rownames(seuobj))
  
  ggfunc <- function(gene,data,x){
    ggplot(data,aes(x = .data[[x]], y = .data[[gene]])) + 
      geom_point(aes(color = cluster), alpha = 1) + 
      geom_smooth(color = "red") + theme_classic() +
      ridgetheme + 
      scale_color_manual(values = colors_list,drop = F) + ggtitle(gene) + 
      theme(legend.position = "none") + xlab("") + ylab("") + 
      guides(colour = guide_legend(override.aes = list(size=10)))
  }
  plot_grid(
    plotlist = lapply(genes, ggfunc, data, x = xaxis, ...),
    labels = labels, ncol = ncol, nrow = ncol, label_size = 18
  )
}

# seurat4 函数
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
# data slot是标准化后的数据 counts slot是没有标准化后的数据
get_data_table <- function(seuobj, highvar = T, type = "counts") {
    if (type == "cpm") {
        countsMat <- as.matrix(GetAssayData(seuobj, slot = "counts"))
    } else {
        countsMat <- as.matrix(GetAssayData(seuobj, slot = type))
    }
    if (highvar == T) {
        # genelist <- GetAssayData(seuobj, slot = "var.features") # 在高版本seurat中弃用
        genelist <- seuobj@assays[["SCT"]]@var.features
        countsMat <- countsMat[genelist, ] # 过滤高变异基因
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

GO_dotplot <- function(gene_list, OrgDb = org.Hs.eg.db, ...) {
    genes <- bitr(gene_list,
        fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"),
        OrgDb = OrgDb
    )
    enrich.go <- enrichGO(
        gene = genes$ENTREZID, # 基因列表文件中的基因名称
        OrgDb = OrgDb,
        keyType = "ENTREZID",
        ont = "ALL", # 可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
    )
    dotplot(enrich.go, title = paste("markers", "GO"), showCategory = 15) + theme_classic() + theme(text = element_text(colour = "black", size = 16), 
                                                                                                  plot.title = element_text(size = 16,color="black",hjust = 0.5),
                                                                                                  axis.title = element_text(size = 16,color ="black"), 
                                                                                                  axis.text = element_text(size= 16,color = "black"))
}

# 3-dim featureplot
surfaceplot <- function(gene, seuobj, x_seq, y_seq, sep = 0.3, contours = F, alpha = 1, baseline = 0, ...) {
    embedding <- FetchData(object = seuobj, vars = c("UMAP_1", "UMAP_2", gene))
    colnames(embedding) <- c("UMAP_1", "UMAP_2", "feature")
    sep <- x_seq[2] - x_seq[1]
    position <- matrix(nrow = length(x_seq), ncol = length(y_seq))

    # umap重新采样，计算z轴位置
    for (q in c(1:length(y_seq)))
    {
        for (p in c(1:length(x_seq)))
        {
            temp <- embedding$feature[embedding$UMAP_1 > x_seq[p] - sep / 2 & embedding$UMAP_1 < x_seq[p] + sep / 2 &
                embedding$UMAP_2 > y_seq[q] - sep / 2 & embedding$UMAP_2 < y_seq[q] + sep / 2]
            position[q, p] <- mean(temp)
        }
    }
    position[is.na(position)] <- baseline

    # 滑动平均
    for (p in c(1:length(x_seq)))
    {
        for (q in c(2:(length(y_seq) - 1))) {
            position[p, q] <- 0.25 * position[p, q - 1] + 0.5 * position[p, q] + 0.25 * position[p, q + 1]
        }
    }

    for (q in c(1:length(y_seq)))
    {
        for (p in c(2:(length(x_seq) - 1))) {
            position[p, q] <- 0.25 * position[p - 1, q] + 0.5 * position[p, q] + 0.25 * position[p + 1, q]
        }
    }

    # 绘图
    if (contours) {
        plot_ly(x = y_seq, y = x_seq, z = t(position), colors = c("lightgrey", "#ff2121"), alpha = alpha) %>%
            add_surface(contours = list(
                z = list(
                    show = TRUE,
                    usecolormap = TRUE,
                    highlightcolor = "#1E90FF",
                    project = list(z = TRUE)
                )
            )) %>%
            layout(scene = list(aspectmode = "manual", aspectratio = list(x = 1, y = 1, z = 0.5)))
    } else {
        plot_ly(x = x_seq + sep / 2, y = y_seq + sep / 2, z = position, colors = c("lightgrey", "#ff212100"), alpha = alpha) %>%
            add_surface() %>%
            layout(scene = list(aspectmode = "manual", aspectratio = list(x = 1, y = 1, z = 0.5)))
    }
}



surfaceplot2 <- function(gene, seuobj, x_seq, y_seq, sep = 0.3, alpha = 1, baseline = 0, z_height=4, ...) {
    embedding <- FetchData(object = seuobj, vars = c("UMAP_1", "UMAP_2", gene,"seurat_clusters"))
    colnames(embedding) <- c("UMAP_1", "UMAP_2", "feature","seurat_clusters")
    
    position <- matrix(nrow = length(x_seq), ncol = length(y_seq))

    # umap重新采样，计算z轴位置
    for (q in c(1:length(y_seq)))
    {
        for (p in c(1:length(x_seq)))
        {
            temp <- embedding$feature[embedding$UMAP_1 > x_seq[p] - sep / 2 & embedding$UMAP_1 < x_seq[p] + sep / 2 &
                embedding$UMAP_2 > y_seq[q] - sep / 2 & embedding$UMAP_2 < y_seq[q] + sep / 2]
            position[q, p] <- mean(temp)
        }
    }
    position[is.na(position)] <- baseline

    # 滑动平均
    for (p in c(1:length(x_seq)))
    {
        for (q in c(2:(length(y_seq) - 1))) {
            position[p, q] <- 0.25 * position[p, q - 1] + 0.5 * position[p, q] + 0.25 * position[p, q + 1]
        }
    }

    for (q in c(1:length(y_seq)))
    {
        for (p in c(2:(length(x_seq) - 1))) {
            position[p, q] <- 0.25 * position[p - 1, q] + 0.5 * position[p, q] + 0.25 * position[p + 1, q]
        }
    }

    # 绘图
    plot_ly(x = x_seq + sep / 2, y = y_seq + sep / 2, z = position, colors = c("#FFFFFF00", "#ff212100"), alpha = alpha) %>%
        add_surface() %>%
        # layout(
        #     scene = list(
        #         aspectmode = "manual",
        #         aspectratio = list(x = 1, y = 1, z = 0.5),
        #         xaxis = list(
        #             title = "UMAP_1",
        #             showgrid = F,
        #             showticklabels=F,
        #             zerolinecolor = "#FFFFFF00"
        #         ), yaxis = list(
        #             title = "UMAP_2",
        #             showgrid = F,
        #             showticklabels=F,
        #             zerolinecolor = "#FFFFFF00"
        #         ), zaxis = list(
        #             title = "expression",
        #             showgrid = F,
        #             zerolinecolor = "#FFFFFF00"
        #         ),
        #         camera = list(eye = list(x = -1.25, y = -1.25, z = 1.25))
        #     )
        # ) %>%
      layout(
        scene = list(
          aspectmode = "manual",
          aspectratio = list(x = 1, y = 1, z = 0.5),
          xaxis = list(
            title = "",
            showgrid = F,
            showticklabels=F,
            zerolinecolor = "#FFFFFF00"
          ), yaxis = list(
            title = "",
            showgrid = F,
            showticklabels=F,
            zerolinecolor = "#FFFFFF00"
          ), zaxis = list(
            title = "", dtick = 1,
            showgrid = F,showticklabels=F,
            zerolinecolor = "#FFFFFF00"
          ),
          camera = list(eye = list(x = -1.25, y = -1.25, z = 1.25))
        )
      ) %>%
        add_trace(
            x = embedding$UMAP_1, y = embedding$UMAP_2, z = z_height, type = "scatter3d",
            inherit = F, split = embedding$seurat_clusters, mode = "markers",
            marker = list(
                symbol = "circle", opacity = 0.2, size = 1.5,
                line = list(color = "#FFFFFF00", width = 1)
            )
        )%>%layout(title = gene, font = list(family = "Arial", size = 25, color = "black"),showlegend = FALSE)
}



# 默认是对于findsubcluster的结果进行
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


## 获得sample info和gene expression
## 分裂小提琴图
func1 <- function(gene, sample, datable){
  data.frame(expr = datable[gene,], sample = sample, gene = gene)
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

## 分裂小提琴图 函数原型
# ggobj <- ggplot(merge_expr,aes(x = gene, y = expr,fill = sample)) +
#   geom_split_violin(trim= F, color="white", scale = "area") + 
#   geom_point(data = Data_summary,aes(x = gene, y= expr), pch=19,
#              position=position_dodge(0.2),size= 1) + #绘制均值位置
#   geom_errorbar(data = Data_summary, aes(ymin = expr-ci, ymax= expr+ci), 
#                 width= 0.05, 
#                 position= position_dodge(0.2), #误差线位置，和均值位置相匹配
#                 color="black",
#                 alpha = 0.7,
#                 size= 0.5) +
#   scale_fill_manual(values = c("#b1d6fb", "#fd9999"))+ 
#   labs(y=("Log2 expression"),x=NULL,title = "Split violin") + 
#   theme_classic()+ mytheme + stat_compare_means(aes(group = sample),
#                                                 label = "p.format",
#                                                 method = "wilcox.test",
#                                                 label.y = max(merge_expr$expr),
#                                                 hide.ns = F)


reg_vector <- function(x){
  paste0("c(",toString(x),")")
}


# for ORA -----------------------------------------------------------------

compute_ratio <- function(x){sapply(x, function(y){parse(text=y) %>% eval()})}
compute_enrichment_factor <- function(x) {return(x@result %>% mutate(enrichment_factor=compute_ratio(GeneRatio)/compute_ratio(BgRatio)) %>%  arrange(desc(enrichment_factor)))}
# 
# GO_signif <- inner_join(df1, df2, by='ID') 
# colnames(GO_signif)
# ggplot(GO_signif %>% filter(ID %in%c("GO:0006487", "GO:0006096", "GO:0006006", "GO:0042593")), aes(x=0, y=factor(Description.x))) + 
#   geom_segment(aes(x=enrichment_factor.y, xend=enrichment_factor.x,  yend=..y..)) + geom_point(aes(x=enrichment_factor.x), color="#4591c7" , size=3) + geom_point(aes(x=enrichment_factor.y), color="#ff5a2c", size=3) + xlim(0,NA) + labs(x="enrichment factor", y="GO term") + theme_bw() + mytheme2
