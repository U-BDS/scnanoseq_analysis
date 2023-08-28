
# This script is a copy of the analysis done for validation of the BLAZE barcode calling tool.
# You, Y., Prawer, Y.D.J., De Paoli-Iseppi, R. et al. Identification of cell barcodes from long-read single-cell RNA-seq with BLAZE. Genome Biol 24, 66 (2023). https://doi.org/10.1186/s13059-023-02907-y

set.seed(1234)

library(ggplot2)
library(Seurat)

#######################
### UMAP DEFINITION ###
#######################
# test
plot_umap <- function(count.matrix, min.features = 200,max.feature = 999999999, npc = 5, cluster_res = 0.7, fig_name = ''){
  
  # init 
  rst_figures <- list()
  rst_table = data.frame()
  
  counts <- read.csv(count.matrix, sep=",")
  
  seurat_object <- CreateSeuratObject(counts = counts, project = "singlecell", min.cells = 3, min.features = 0.5)
  rst_table <- rbind(rst_table, data.frame("Cells"=dim(seurat_object)[2],
                                           "median feature per Cell"=median(seurat_object$nFeature_RNA), row.names = paste0('min features > 0'),check.names = FALSE))
  
  #remove unwanted cells. below are default settings but you can modify these
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > min.features & nFeature_RNA < max.feature) 
  rst_table <- rbind(rst_table, data.frame("Cells"=dim(seurat_object)[2],
                                           "median feature per Cell"=median(seurat_object$nFeature_RNA), row.names = paste0('min features > 0', min.features),check.names = FALSE))
  
  #now you have removed unwanted cells, it is time to normalize the data. By default, Seurat employs a global-scaling normalization method "LogNormalize" that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  #you can alternatively input seurat_object <- NormalizeData(seurat_object) instead of above.
  #we now identify highly variable features in order to determine a subset of features that exhibit high cell-to-cell 
  
  #variation in the dataset.
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  
  #now we apply a linear transformation (scaling) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  #we can visualise both cells and features that define the PCA
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  #this can be visualised using an Elbow Plot
  #rst_figures <- append(rst_figures, ElbowPlot(seurat_object))
  
  #to cluster cells, I have used the number of sig. PCs that I observed in the above plots. The findneighbors function is for constuction of a KNN graph based on euclidean distance in PCA space and refines edge weights between any two cells based on the shared overlap in their local neighborhoods (jaccard similarity). It uses the input of previously defined dimensionality of the dataset.
  seurat_object <- FindNeighbors(seurat_object, dims = 1:npc)
  #now to actually cluster the cells, we apply modularity optimisation techniques (default is Louvain algorithm). The findclusters function contains a resolution parameter which sets the granularity of downstream clustering. Settings are recommended between 0.4-1.2, but may need to be increased for larger datasets.
  seurat_object <- FindClusters(seurat_object, resolution = cluster_res)
  
  #to View metadata
  #seurat_object@meta.data
  #run non-linear dimensional reduction (UMAP/tSNE)
  seurat_object <- RunUMAP(seurat_object, dims = 1:npc)
  
  rst_figures <- append(rst_figures, list(DimPlot(seurat_object, reduction = "umap") + ggplot2::labs(color = "cluster \n(from PCA)", title = '') + ggplot2::theme(text = ggplot2::element_text(size = 10))  )) 
  
  
  rst_figures <- append(rst_figures, list(
    FeaturePlot(seurat_object, reduction = "umap", features = 'nCount_RNA')+ggplot2::labs(color = "UMI count",title = '')+ ggplot2::theme(text = ggplot2::element_text(size = 10)),
    FeaturePlot(seurat_object, reduction = "umap", features = 'nFeature_RNA')+ggplot2::labs(color = stringr::str_wrap("Feature count (isoform/gene)",15),title = '')+ ggplot2::theme(text = ggplot2::element_text(size = 10))
  ))

  plot_pc <- ElbowPlot(seurat_object)+ggplot2::labs(title = 'SD explained by each PC') + ggplot2::theme(text = ggplot2::element_text(size = 10))
  plot_umap <- gridExtra::grid.arrange(plot_pc, gridExtra::tableGrob(rst_table), rst_figures[[1]], rst_figures[[2]], rst_figures[[3]], ncol=2, top=grid::textGrob(fig_name))
  list(plot_umap, seurat_object)
  
  seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, logfc.threshold = 0.25)
  # Group All Markers by cluster id
  grouped_markers <- seurat_object.markers %>%
    group_by(cluster)
  
  
}


##################
### INPUT DATA ###
##################

# BLAZE nonQ20
transcript_mtx <- "input/blaze/trans_count_B+F_non_q20_rebase.csv"
gene_mtx <- "input/blaze/gene_count_B+F_non_q20_rebase.csv"


################
### ANALYSES ###
################

plots <- plot_umap(transcript_mtx,
                   min.features = 200,
                   max.feature = 999999999,
                   npc = 5,
                   cluster_res = 0.7,
                   fig_name = 'BLAZE+FLAME(Transcript counts, non-Q20 Gridion data)')


plots