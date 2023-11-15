# bulkRNA-seq wrappers ----------------------------------------------------------
# create: 5th Nov., 2021
# last updated: 2023/7/21

library(dplyr)
suppressMessages(library(DESeq2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(Mfuzz))
read_matrix<-function(file_name) ##from the given file import count matrix,return data matrix
{
  data<-read.csv(file_name, check.names = FALSE)
  datatrix<-as.matrix(data[,c(-1)])
  rownames(datatrix)<-data$Geneid
  #head(datatrix) #show first 6
  return(datatrix)
}

make_grouping_file<-function(datatrix, Condition=NULL, levels=NULL, dec='-')
{
  if (!is.na(levels)){
    # if levels given, do the automatic recognition
    temp <- strsplit(colnames(datatrix), dec, fixed=TRUE)
    Condition <- sapply(temp, function(x){x[x %in% levels]})
  }
  grouping_file<-data.frame(row.names <- colnames(datatrix),Condition)
  return(grouping_file)
}
prefilter<-function(datatrix)
{
  datatrix <- datatrix[rowSums(datatrix)>= 10,]
  datatrix <- datatrix[rowSums(datatrix<1)<ncol(datatrix)/3, ] # update 0721: filter out genes with limited expression
  return(datatrix)
}

make_dds<-function(dataset,grouping_file)
{
  dmatrix<-DESeqDataSetFromMatrix(countData=dataset,colData=grouping_file,
                                  design=~Condition)
  return(dmatrix)
}

suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))
quality_control<-function(analyseddata)
{
  vsdtrans <<- vst(analyseddata, blind=FALSE)
  #rldtrans <- rlog(analyseddata, blind=FALSE)
  ntdtrans <<- normTransform(analyseddata)
  sampleDists <- dist(t(assay(vsdtrans)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsdtrans$row.names....colnames.datatrix.)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
  ggobj<-plotPCA(vsdtrans, intgroup=c("Condition"),returnData=FALSE)
  plot(ggobj+ theme_bw() + geom_text_repel(label = ggobj[["data"]][["name"]]))
}

DEGs <- function(file_name,FC=2,padj=0.001)#return differently expressed genes
{
  sf <- as.data.frame(read.csv(file_name))
  sf <- mutate_all(sf, ~replace(., is.na(.), 0)) # (important) update 230721: replace NA in pvalue and p.adj
  sf$threshold = factor(ifelse(sf$padj < padj & abs(sf$log2FoldChange) >= FC, 
                               ifelse(sf$log2FoldChange>= FC ,'Up','Down'),'NoSignifi'),
                        levels=c('Up','Down','NoSignifi'))
  
  orderedsf<- sf %>% arrange(desc(log2FoldChange))
  upgene<-orderedsf$X[which(orderedsf$threshold=="Up")]
  downgene<-orderedsf$X[which(orderedsf$threshold=="Down")]
  res = list(upgene,downgene)
  return(res)
}

vo_plot <- function(file_name, FC = 3, padj = 0.001) {
  sf <- as.data.frame(read.csv(file_name))
  sf <- NaRV.omit(sf)
  sf$threshold <- factor(ifelse(sf$padj < padj & abs(sf$log2FoldChange) >= FC,
                                ifelse(sf$log2FoldChange >= FC, "Up", "Down"), "NoSignifi"
  ),
  levels = c("Up", "Down", "NoSignifi")
  )
  
  orderedsf <- sf[order(sf$log2FoldChange), ]
  
  uptopsf <<- head(sf$X[which(sf$threshold == "Up")], 10)
  downtopsf <<- head(sf$X[which(sf$threshold == "Down")], 10)
  
  topssf <- c(as.character(uptopsf), as.character(downtopsf))
  
  
  orderedsf$Label <- ""
  orderedsf$Label[match(topssf, orderedsf$X)] <- topssf
  
  ggplot(orderedsf, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
    geom_point() +
    scale_color_manual(values = c("#CC0000", "#BBBBBB", "#2f5688")) +
    geom_text_repel(
      data = orderedsf[orderedsf$padj < padj & abs(orderedsf$log2FoldChange) >= FC, ],
      aes(label = Label),
      size = 3,
      segment.color = "black", show.legend = FALSE
    ) +
    theme_bw() +
    theme(
      legend.title = element_blank()
    ) +
    ylab("-log10 (p-adj)") +
    xlab("log2 (FoldChange)") +
    geom_vline(xintercept = c(-FC, FC), lty = 3, col = "black", lwd = 0.5) +
    geom_hline(yintercept = c(0, 3), lty = 3, col = "black", lwd = 0.5)
}

heat_map <- function() {
  selecteddata <- assay(ntdtrans)[c(uptopsf, downtopsf), ]
  pheatmap(selecteddata,
           display_numbers = FALSE, number_color =
             "black", cluster_rows = FALSE, cluster_cols = FALSE, scale = "row"
  )
}