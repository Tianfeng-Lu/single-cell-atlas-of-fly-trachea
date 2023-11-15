scVolcano <- function(diffData, topGeneN=5, orderBy="avg_log2FC", plotTitle = NULL, 
                      backH=0.1, annoH = 0.5, myMarkers=NULL, clusterOrder=NULL, seed=143, cluster_colors=colors_list[-1]){
  set.seed(seed)
  diff.marker <- diffData %>% mutate(
    pattern = ifelse(avg_log2FC>0, "Up", "Down")
  )
  
  if (!is.null(clusterOrder)) {
    diff.marker$cluster <- factor(diff.marker$cluster, levels = clusterOrder)
  }
  
  back.data <- purrr::map_df(unique(diff.marker$cluster), function(x){
    tmp <- diff.marker %>% filter(cluster == x)
    new.tmp <- data.frame(
      cluster = x,
      min = min(tmp$avg_log2FC) - backH,
      max = max(tmp$avg_log2FC) + backH
    )
    return(new.tmp)
  })
  
  top.marker.tmp <- diff.marker %>% group_by(cluster)
  top.marker.max <- top.marker.tmp %>% 
    slice_max(n = topGeneN, order_by = get(orderBy))
  top.marker.min <- top.marker.tmp %>%
    slice_min(n = topGeneN, order_by = get(orderBy))
  top.marker <- rbind(top.marker.max, top.marker.min)
  
  if (!is.null(myMarkers)) {
    top.marker <- diff.marker %>%
      filter(gene %in% myMarkers)
  }else{
    top.marker <- top.marker
  }
  
  diff.marker$topN <- ifelse(
    paste0(diff.marker$cluster, diff.marker$gene) %in% paste0(top.marker$cluster, top.marker$gene), 
    diff.marker$gene, "")
  
  diff.marker$topNShape <- ifelse(
    paste0(diff.marker$cluster, diff.marker$gene) %in% paste0(top.marker$cluster, top.marker$gene), 
    1, "")
  
  diff.marker$topNCol <- ifelse(
    paste0(diff.marker$cluster, diff.marker$gene) %in% paste0(top.marker$cluster, top.marker$gene), 
    "black", "white")
  
  diff.marker$topNAlpha <- ifelse(
    paste0(diff.marker$cluster, diff.marker$gene) %in% paste0(top.marker$cluster, top.marker$gene), 
    1, 0)
  

  
  diff.marker.xj <- diff.marker
  diff.marker.xj$xj <- jitter(as.numeric(factor(diff.marker$cluster)), amount = 0.4)
  
  p1 <- ggplot(data = diff.marker.xj, aes(x = xj, y = avg_log2FC)) +
    geom_col(data = back.data, aes(x = cluster, y = min), fill = "grey93") +
    geom_col(data = back.data, aes(x = cluster, y = max), fill = "grey93") +
    geom_point(aes(color = pattern)) +
    geom_point(
      data = diff.marker.xj, 
      aes(x = xj, y = avg_log2FC, shape = topNShape), 
      size = 1.8, stroke = 0.8,
      color = diff.marker.xj$topNCol,
      alpha = diff.marker.xj$topNAlpha,
      show.legend = FALSE
    ) +
    scale_shape_manual(
      values = c(0, 1)
    ) +
    xlab(label = "Clusters") +
    ylab(label = "log2FoldChange") +
    scale_color_manual(
      values = c(Up = "#fd9999", Down = "#b1d6fb"), limits=c("Up", "Down")
    ) +
    scale_y_continuous(
      n.breaks = 6
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.background = element_blank()
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 5)
      )
    )
  # p1 
  
  p2 <- p1 + geom_tile(
    aes(x = cluster, y = 0, fill = cluster),
    # color of tile edge
    color = "black",
    # height of tiles
    height = annoH,
    alpha = 0.5, 
    show.legend = FALSE
  ) +
    ggrepel::geom_text_repel(
      data = diff.marker.xj,
      aes(x = xj, y = avg_log2FC, label = topN), min.segment.length=0.5,
      max.overlaps = 1000, force = 1, force_pull = 2, max.iter=1e4, size=4
    ) + scale_fill_manual(values = cluster_colors)
  # p2
  
  p3 <- p2 + 
    ggtitle(label = "") +
    geom_text(aes(x = ., y = 0, label = .), data=(unique(diffData$cluster) %>% data.frame()), size = 4) +
    theme_bw()+
    theme(
      axis.line.x = element_blank(), 
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(color='black', size=10)
      
    )
  
  return(p3) 
}
