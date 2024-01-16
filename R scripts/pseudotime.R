# pseudotime analysis

library(monocle3)
source("./cleanRwrappers.R")

ctrl_filtered <- readRDS("./ctrl_filtered.rds")
PC_DB_DT <- subset(ctrl_filtered, idents=c("PC","DB","DT"))

cds <- as.cell_data_set(PC_DB_DT)
cds <-  reduce_dimension(cds)
cds@int_colData@listData[["reducedDims"]]@listData[["PCA"]] <- PC_DB_DT@reductions[["pca"]]@cell.embeddings
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- PC_DB_DT@reductions[["umap"]]@cell.embeddings
rowData(cds)$gene_short_name <- row.names(rowData(cds))
colData(cds)@listData[["seurat_clusters"]] <- PC_DB_DT@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- PC_DB_DT@active.ident
cds@clusters@listData[["UMAP"]][["partitions"]] <-  PC_DB_DT@active.ident


# run Monocle3
cds <- cluster_cells(cds,resolution = 1e-5)
cds <- learn_graph(cds,use_partition = F)
cds <- choose_cells(cds)
cds <- order_cells(cds)


# identify genes differentially regulated in pseudotime
pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=48)
pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=0.001)
names(gene_module_df)[1] <- "gene_short_name"
pr_merge <- merge(pr_test_res, gene_module_df, by = "gene_short_name")


top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


pseudotime_table <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime_table <- pseudotime_table[rownames(PC_DB_DT@meta.data)] %>% data.frame()

#add feature to seurat object: pseudotime
PC_DB_DT$pseudotime <- pseudotime_table
# saveRDS(PC_DB_DT, "./PC_DB_DT_seurat.rds")

### end of Pseudotime ananlysis----

