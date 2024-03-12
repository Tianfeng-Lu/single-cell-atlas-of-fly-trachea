# A single-cell atlas of fly trachea
Analysis and visualization code for the following article:
Li Y, Lu T, Dong P, Chen J, Zhao Q, Wang Y, Xiao T, Wu H, Zhao Q and Huang H. A single-cell atlas of *Drosophila* trachea reveals glycosylation-mediated Notch signaling in cell fate specification. *Nature Communications* **15**, 2019 (2024). [https://doi.org/10.1038/s41467-024-46455-w](https://doi.org/10.1038/s41467-024-46455-w)

The processed data files are available with [doi: 10.5281/zenodo.10139562](https://doi.org/10.5281/zenodo.10139562) for reproducing main figures in the article and further exploration. You can start anywhere after loading `.rds` files with the same name as dependent variables.

## Source data
### download from zenodo
`download_data.r` provides a pipeline for downloading the raw data and processed data from zenodo.


## Single cell RNA-seq analysis
### Preprocess
`preprocess.R` transforms CellRanger output to Seurat objects with cell type annotation. It keeps main components in trachea and split to filtered ctrl and HSD data which will be used in the downstream analyses.  

Saved files: 
1. fly_merge_filtered.rds 
2. ctrl_filtered.rds
3. HSD_filtered.rds

```Note: Clustering and visualization results might be slightly different.```

### Overview: Gene expression and cell population visualization 
`10x figures.r` provides an overview for our scRNA-seq data, in which we explore gene expression heterogeneity in different cell population.
It needs ctrl_filtered and HSD_filtered in the R environment.


### Analysis for differential expression genes
DEG analysis is implemented in `DEG analysis.R`.
figures:  
    single cell volcano plot compares gene expression in  within each cell cluster   
    ORA shows Notch related pathways are enriched in ND vs HSD

### Pseudotime
`pseudotime.R` performs pseudotime analysis for PC_DB_DT branches via monocle3 software. `pseudotime_visualization.Rmd` is the visualization code in which gene expression~pseudotime and pseudotime heatmap are plotted.


Note: Both `pseudotime.R` and `pseudotime_visualization.Rmd` require `monocle3` to be installed. You might also need to interact with python environment via `reticulate` R package.
```
devtools::install_github("cole-trapnell-lab/monocle3") # install monocle3 from source
```

Saved files:
1. PC_DB_DT_seurat.rds
2. PC_DB_DT_cds.rds

### Subclustering
To further explore gene expression heterogeneity within a cell population, we exploit subclustering to clusters of interest in `subclustering.r`. 

### RNA velocity
RNA velocity analysis via python package `scvelo`, input file is `ctrl_velo.h5ad`.
The code is available in `scvelo.ipynb`.

## Bulk RNA-seq analysis
`./bulk RNA-seq/CodeForTimeDifferentialAnalysis.Rmd` provides a pipeline to process bulk RNA-seq raw data via DESeq2. And then performs ORA for DEGs upregulated when L3 -> 0hr. 


## Utils
### UpSetR
Upset plot in R with modification to plot percentage rather than counts.
To enable it
```
devtools::load_all("./UpSetR-master", export_all = TRUE)
```

### scVolcano
Volcano plot for each cluster in single cell data.

### cleanRwrappers
`cleanRwrappers.R` contains some useful functions for visualization and coloring.

### bulkRNAseqWrappers
`./bulk RNA-seq/bulkRNAseqWrappers.R` contains useful functions for bulk RNAseq data processing and visualization.
