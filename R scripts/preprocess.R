# preprocess: clustering and annotation
# preprocessed data file: ctrl_filtered.rds, HSD_filtered.rds

source("./cleanRwrappers.R")

# load raw data
fly_merge <- CreateSeuratObject(Read10X("./fly_merge/"), names.field = 2, names.delim = "-", project = "fly_merge", min.cells = 10, min.features = 200)

fly_merge <- PercentageFeatureSet(fly_merge, pattern = "^mt:", col.name = "percent.mt")

VlnPlot(fly_merge,c("percent.mt"))/
  VlnPlot(fly_merge,c("nFeature_RNA"))/
  VlnPlot(fly_merge,c("nCount_RNA"))

fly_merge <- subset(fly_merge, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & 
                      nCount_RNA > 1000 &  nCount_RNA < 20000 & percent.mt < 10)

fly_merge <- fly_merge %>% SCTransform(vars.to.regress = "percent.mt", verbose = F) %>% 
  RunPCA() %>% FindNeighbors(dims = 1:20) %>% 
  RunUMAP(dims = 1:20) %>% 
  FindClusters(resolution = 0.1)


# cell type annotation
umapplot(fly_merge)
Idents(fly_merge) <- fly_merge$seurat_clusters
levels(Idents(fly_merge)) <- c("DB","DT","SB","PC","Fat Body","Muscle","TC","VB","LT","NE","GB")

selected_cells <- CellSelector(plot = DimPlot(fly_merge, reduction = "umap"))
fly_merge <- SetIdent(fly_merge, cells = selected_cells, "ASP")
fly_merge$celltype <- Idents(fly_merge)

# only keep main components in trachea
fly_merge_filtered <- subset(fly_merge, 
                             idents = c("DB","DT","SB","PC","TC","VB","LT","GB","ASP"))
fly_merge_filtered$celltype <- droplevels(fly_merge_filtered$celltype) 

fly_merge_filtered <- fly_merge_filtered %>% 
  RunPCA() %>% FindNeighbors(dims = 1:20) %>% 
  RunUMAP(dims = 1:20)
umapplot(fly_merge_filtered,split.by = "orig.ident")

fly_merge_filtered <- NormalizeData(fly_merge_filtered, assay = "RNA") %>% ScaleData(assay = "RNA")

# split ctrl-HSD
ctrl_filtered <- subset(fly_merge_filtered, orig.ident == 1)
HSD_filtered <- subset(fly_merge_filtered, orig.ident == 2 | orig.ident == 3)


