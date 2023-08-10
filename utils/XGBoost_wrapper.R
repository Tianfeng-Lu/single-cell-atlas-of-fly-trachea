library(xgboost)
library(Matrix)
library(mclust)
library(tidyverse)
library(SHAPforxgboost)
library(lambda.r)
XGBoost_train_from_seuobj <- function(seuobj, is_highvar = T, test_ratio = 0.3, seed = 7)
{ 
  ## set test_ratio to 0 to avoid extracting test from dataset
  set.seed(seed)
  seuobj_label <- as.numeric(as.character(Idents(seuobj)))
  if(is.na(seuobj_label[1])) # check vaild Idents
  {
    warning("Please ensure that seurat idents are in numeric forms")
    seuobj_label <- as.numeric(Idents(seuobj))-1 # transform character idents to numeric
  }
  # colnames(seuobj_data) <- NULL
  seuobj_data <- get_data_table(seuobj, highvar = T, type = "data")
  xgb_param <- list(eta = 0.2, max_depth = 6, 
                    subsample = 0.6,  num_class = length(table(Idents(seuobj))),
                    objective = "multi:softprob", eval_metric = 'mlogloss')
  
  if(test_ratio == 0) {
    seuobj_train_data <- list(data = t(as(seuobj_data,"dgCMatrix")), label = seuobj_label) 
    # use whole dataset as train data
    seuobj_train <- xgb.DMatrix(data = seuobj_train_data$data,label = seuobj_train_data$label)
    bst_model <- xgb.train(xgb_param, seuobj_train, nrounds = 100, verbose = 0)
  } else {
    index <- c(1:dim(seuobj_data)[2]) %>% sample(ceiling(test_ratio*dim(seuobj_data)[2]), replace = F, prob = NULL)
    seuobj_train_data <- list(data = t(as(seuobj_data[,-index],"dgCMatrix")), label = seuobj_label[-index])
    seuobj_test_data <- list(data = t(as(seuobj_data[,index],"dgCMatrix")), label = seuobj_label[index])
    seuobj_test <- xgb.DMatrix(data = seuobj_test_data$data,label = seuobj_test_data$label)
    seuobj_train <- xgb.DMatrix(data = seuobj_train_data$data,label = seuobj_train_data$label)
    watchlist <- list(train = seuobj_train, eval = seuobj_test)
    bst_model <- xgb.train(xgb_param, seuobj_train, nrounds = 100, watchlist, verbose = 0)
  }
  return(bst_model)
}

show_train_loss <- function(bst_model, nrounds = 100) #when $ test_ratio \neq 0 $ show loss in watchlist
{
  eval_loss <- bst_model[["evaluation_log"]][["eval_mlogloss"]]
  plot_ly(data.frame(eval_loss), x = c(1:nrounds), y = eval_loss) %>% 
    add_trace(type = "scatter", mode = "markers+lines", 
              marker = list(color = "black", line = list(color = "#1E90FFC7", width = 1)),
              line = list(color = "#1E90FF80", width = 2)) %>% 
    layout(xaxis = list(title = "epoch"),yaxis = list(title = "eval_mlogloss"), 
           title = "train_loss", font = list(family = "Arial", size = 25, color = "black"))
}

XGBoost_predict_from_seuobj <- function(seuobj, bst_model, is_highvar = T, seed = 7, celltype_assign = 2, return_confuse_matrix = F)
{
  #return a updated seurat object with new metadata named confidence and projected_idents
  seuobj_label <- as.numeric(as.character(Idents(seuobj)))
  if(!is.null(which(is.na(seuobj_label)))) # check vaild Idents
  {
    warning("Please ensure that seurat idents are in numeric forms")
    seuobj_label <- as.numeric(Idents(seuobj))-1 
  }
  temp <- get_data_table(seuobj, highvar = T, type = "data")
  seuobj_data <- matrix(data = 0, nrow = bst_model$nfeatures, ncol = length(colnames(temp)), 
                        byrow = FALSE, dimnames = list(bst_model[["feature_names"]],colnames(temp)))
  intersect_features <- intersect(bst_model[["feature_names"]], rownames(temp))
  seuobj_data[intersect_features,] <- temp[intersect_features,]
  rm(temp)
  
  
  # colnames(seuobj_data) <- NULL
  seuobj_test_data <- list(data = t(as(seuobj_data,"dgCMatrix")), label = seuobj_label)
  seuobj_test <- xgb.DMatrix(data = seuobj_test_data$data,label = seuobj_test_data$label)
  
  #预测结果
  predict_seuobj_test <- predict(bst_model, newdata = seuobj_test)
  predict_prop_seuobj <<- matrix(data=predict_seuobj_test, nrow = bst_model[["params"]][["num_class"]], 
                                 ncol = dim(seuobj_test_data$data)[1], byrow = FALSE, 
                                 dimnames = list(as.character(0:(bst_model[["params"]][["num_class"]]-1)),
                                                 colnames(seuobj)))
  
  # predict cell types
  if(celltype_assign == 1){
    seuobj_res <- apply(predict_prop_seuobj,2,ident_assignfunc,rownames(predict_prop_seuobj))
  }else if(celltype_assign == 2){
    seuobj_res <- apply(predict_prop_seuobj,2,ident_assignfunc2,rownames(predict_prop_seuobj))
  }else if(celltype_assign == 0){
    seuobj_res <- apply(predict_prop_seuobj,2,ident_assignfunc0,rownames(predict_prop_seuobj))
  }else{
    stop("celltype_assign invalid")
  }
  
  confuse_matrix <- table(seuobj_test_data$label, seuobj_res, dnn=c("true","pre"))
  
  print(paste('ARI =',adjustedRandIndex(seuobj_res, seuobj_test_data$label)))
  
  seuobj <- AddMetaData(seuobj, data.frame(t(predict_prop_seuobj), stringsAsFactors=F))
  
  #save and update seurat object
  seuobj$projected_idents <- factor(seuobj_res, levels =
                                      c(as.character(0:(bst_model[["params"]][["num_class"]]-1)),"unassigned"))

  if(return_confuse_matrix){
    return(confuse_matrix)
  }else{
    print("return a seurat object with meta.data'X1'~'Xn'")
    print("return a seurat object with meta.data'projected_idents'")
  return(seuobj)
  }
}

## assign cell type predicted via tree models, consider confidence
ident_assignfunc <- function(s, ident) {
  if (max(s) > 1.5 / length(ident)) {
    return(ident[which(s == max(s))])
  } else {
    return("unassigned")
  }
}

ident_assignfunc2 <- function(s, ident)
  # confidence : max - 2th_max > 0.4
{
  if (max(s) - max(s[s!=max(s)]) > 0.4) {
    return(ident[which(s == max(s))])
  } else {
    return("unassigned")
  }
}

ident_assignfunc0 <- function(s, ident)
  # softmax
{
    return(ident[which(s == max(s))])
}

project2ref_celltype <- function(query_seuobj, ref_seuobj, 
                                 query_labels = "projected_idents",
                                 ref_labels = c("seurat_clusters","Classification1")) 
  # transfer lables: add ref_celltype to meta.data in query seurat object
  # ref_labels assign the mapping between numeric idents and celltype idents
{
  identmap <- 0:(length(levels(ref_seuobj@meta.data[[ref_labels[1]]]))-1)
  ## build mapping between numeric labels and ref labels
  names(identmap) <- levels(ref_seuobj@meta.data[[ref_labels[2]]]) 
  
  df <- query_seuobj@meta.data[[query_labels]]
  lambda(x,identmap) %:=% ifelse(x=="unassigned",x,names(identmap[identmap == x]))
  levels(df) <- lapply(levels(df), lambda, identmap) %>% as.character() # permute cell labels via `identmap`
  query_seuobj$ref_celltype <- df
  return(query_seuobj)
}

## suitable for current ident to train bst_model
project2ref_celltype2 <- function(query_seuobj, ref_seuobj, 
                                 query_label = "projected_idents") 
{
  df <- Idents(query_seuobj) 
  Idents(query_seuobj) <- query_seuobj@meta.data[[query_label]]
  t1 <- levels(Idents(query_seuobj))
  t2 <- levels(Idents(ref_seuobj))
  if(length(t1) == length(t2)){ # ID convention
  levels(Idents(query_seuobj)) <- c(levels(Idents(ref_seuobj)))
  }else{
    levels(Idents(query_seuobj)) <- c(levels(Idents(ref_seuobj)),"Unassigned")
  }
  query_seuobj$ref_celltype <- Idents(query_seuobj)
  Idents(query_seuobj) <- df
  return(query_seuobj)
}

library(SingleCellExperiment)
library(scmap)

mkref_scmap_from_seuobj <- function(ref_seuobj){
  # return a ref sce object
  ref_sce <- as.SingleCellExperiment(ref_seuobj)
  logcounts(ref_sce) <- log2(counts(ref_sce) + 1)
  
  counts(ref_sce) <- as.matrix(counts(ref_sce))
  logcounts(ref_sce) <- as.matrix(logcounts(ref_sce))
  
  rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
  ref_sce <- ref_sce[!duplicated(rownames(ref_sce)), ]
  ref_sce <- selectFeatures(ref_sce, suppress_plot = T) %>% indexCell()
  
  return(ref_sce)
}

query_scmap_from_refsce <- function(query_seuobj, ref_sce, ref_labels = 'Classification1'){
  #return a updated seurat object with new metadata named scmap_idents判断对象是否存在
  
  query_sce <- as.SingleCellExperiment(query_seuobj)
  logcounts(query_sce) <- log2(counts(query_sce) + 1)
  
  counts(query_sce) <- as.matrix(counts(query_sce))
  logcounts(query_sce) <- as.matrix(logcounts(query_sce))
  
  rowData(query_sce)$feature_symbol <- rownames(query_sce)
  query_sce <- query_sce[!duplicated(rownames(query_sce)), ]
  query_sce <- selectFeatures(query_sce, suppress_plot = T) %>% indexCell()
  
  scmapCell_results <- scmapCell(query_sce, list(ref = metadata(ref_sce)$scmap_cell_index))
  scmapCell_clusters <- scmapCell2Cluster(scmapCell_results,
                                          list(as.character(colData(ref_sce)[[ref_labels]])))
  
  query_seuobj$scmap_idents <- data.frame(scmapCell_clusters$scmap_cluster_labs[,"ref"],
                                          row.names = colnames(query_seuobj))
  return(query_seuobj)
}

upsetplot_celltype <- function(seuobj, cell_type){
  # draw upset plot for the given cell type
  df <- cbind(seuobj[["Classification1"]],seuobj[["ref_celltype"]],seuobj[["scmap_idents"]])
  df <- df[df$Classification1 == cell_type | df$ref_celltype == cell_type | df$scmap_idents == cell_type,]
  li <- table(df) %>% as.data.frame() #获得含有SMC1的frequency
  li <- li[!(li$Freq < 5),] #删除frequency<5的行
  li <- li[order(li$Freq,decreasing = T),]
  dd <- data.frame(li,index = as.character(1:nrow(li)))
  
  dd$index <- factor(dd$index,levels = 1:length(levels(dd$index)))
  dd$Freq <- NULL
  dd <- reshape2::melt(dd,id.var = "index")
  colnames(dd) <- c("index","type","name")
  
  dt <- data.frame(Freq = li$Freq,index = as.character(1:nrow(li)))
  dt$index <- factor(dt$index,levels = 1:length(levels(dt$index)))
  dt$col <- dt$Freq>(sum(dt$Freq)/10) #set the color to red for those frequency > 0.1
  
  p1 <- ggplot(dd)+geom_point(mapping = aes(x = index, y = type, color = name), size = 6) + mytheme2 + scale_color_manual(values = colors_list, drop = F) + theme(axis.ticks.x = element_blank(),axis.title.x = element_blank())
  p2 <- ggplot(dt,aes(x = index, y = Freq, fill = col))+geom_bar(stat = "identity") + mytheme2 + theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(), legend.position = 'none') + scale_fill_manual(values = c("black","red"))
  
  plot <- cowplot::plot_grid(p2,p1,ncol = 1,align = 'v',rel_heights = c(2,1))
}


## unsupervised version
