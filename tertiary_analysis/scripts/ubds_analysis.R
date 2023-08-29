set.seed(1234)

# Set the current directory to the dir above the current directory to access inputs easier
#curr_dir <- getwd()
#setwd(paste0(curr_dir,"/.."))

library(ggplot2)
library(Seurat)
library(dplyr)
library(grid)
library(sctransform)
library(biomaRt)
library(clustree)

# Install the UBDS R scripts
lapply(list.files("./scripts/R"), FUN = function(x) source(paste0("./scripts/R/", x)))

#### CONSTANT DEFS ####

feature_type <- 'transcript'

# The path to the input file
in_file <- paste0("input/isoquant/q20.", feature_type, "_counts.tsv")

# The name of the project to use in seurat
project_name <- "isoquant_q20_seurat"

# The prefix to output files
output_prefix <- paste0("./output/isoquant/q20/", feature_type, "/")

#### PREPARE INPUT DATA ####

dir.create(file.path(output_prefix), recursive = TRUE)

# Read in the matrices
cell_bc_matrix <- read.table(in_file,
                             sep="\t",
                             header = TRUE,
                             row.names = 1)

# Create the seurat object using RNA Assay
seurat_obj <- CreateSeuratObject(counts = cell_bc_matrix,
                                 assay = "RNA",
                                 min.cells = 1,
                                 min.features = 1,
                                 project = project_name)

#### GENERATE QC'S ####

# MITOCHONDRIAL PERCENT
seurat_obj <- PercentageFeatureSet(object = seurat_obj,
                                   pattern = "^mt-",
                                   col.name = "percent.mt",
                                   assay = "RNA")

# VIOLIN PLOTS
pre_filter_violin <- paste0(output_prefix, "violin_prefilter.pdf")

pdf(file = pre_filter_violin,
    width = 8,
    height = 6)

plot(
  VlnPlot(object = seurat_obj,
          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
          ncol = 3,
          assay = "RNA")
  )

dev.off()

# FEATURE SCATTER PLOTS
pre_filter_scatter <- paste0(output_prefix, "scatter_prefilter.pdf")

pdf(file = pre_filter_scatter,
    width = 8,
    height = 6)

plot1 <- FeatureScatter(object = seurat_obj,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.mt")

plot2 <- FeatureScatter(object = seurat_obj,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")

plot(plot1 + plot2) 

dev.off()

# DENSITY PLOTS - nFeature_RNA unscaled
pre_filter_feature_dens <- paste0(output_prefix, "nFeature_unscaled_prefiltered_dens.pdf")

pdf(file = pre_filter_feature_dens,
    width = 8,
    height = 6)

dens_plot_nfeature <- plotSingleCellDensity(input_obj = seurat_obj,
                      metadata_variable = "nFeature_RNA",
                      group.by = "orig.ident",
                      scale_x_axis = FALSE) +
  geom_vline(xintercept = 200)

plot(dens_plot_nfeature)
dev.off()

# DENSITY PLOTS - nCount_RNA unscaled
pre_filter_feature_dens <- paste0(output_prefix, "nCount_unscaled_prefiltered_dens.pdf")

pdf(file = pre_filter_feature_dens,
    width = 8,
    height = 6)

dens_plot_ncount <- plotSingleCellDensity(input_obj = seurat_obj,
                      metadata_variable = "nCount_RNA",
                      group.by = "orig.ident",
                      scale_x_axis = FALSE)  +
  geom_vline(xintercept = 500)

plot(dens_plot_ncount)
dev.off()

# DENSITY PLOTS - nFeature_RNA scaled
pre_filter_feature_scaled_dens <- paste0(output_prefix, "nFeature_scaled_prefiltered_dens.pdf")

pdf(file = pre_filter_feature_scaled_dens,
    width = 8,
    height = 6)

dens_plot_nfeature <- plotSingleCellDensity(input_obj = seurat_obj,
                                            metadata_variable = "nFeature_RNA",
                                            group.by = "orig.ident",
                                            scale_x_axis = TRUE) +
  geom_vline(xintercept = 200)

plot(dens_plot_nfeature)
dev.off()

# DENSITY PLOTS - nCount_RNA scaled
pre_filter_feature_dens_scale <- paste0(output_prefix, "nCount_scaled_prefiltered_dens.pdf")

pdf(file = pre_filter_feature_dens_scale,
    width = 8,
    height = 6)

dens_plot_ncount <- plotSingleCellDensity(input_obj = seurat_obj,
                                          metadata_variable = "nCount_RNA",
                                          group.by = "orig.ident",
                                          scale_x_axis = TRUE)  +
  geom_vline(xintercept = 500)

plot(dens_plot_ncount)
dev.off()

# Scatter Plot All 3 Variables

mt_prefilter <- FetchData(object = seurat_obj,
                          vars = c("orig.ident", 
                                   "nCount_RNA", 
                                   "nFeature_RNA", 
                                   "percent.mt")) %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point() + 
  scale_colour_gradient(low = "grey90", high = "black", limits = c(0,100)) +
  stat_smooth(method = lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident, scales = "free") +
  theme_classic()

mt_prefilter_name <- paste0(output_prefix, "mt_prefilter.pdf")
pdf(mt_prefilter_name,
    width = 8,
    height = 6)
plot(mt_prefilter)
dev.off()

# Save cell count prior to filtering
cell_count <- FetchData(object = seurat_obj, vars = "orig.ident") %>%
              tally(name = "Cells_before_filtering")

write.csv(cell_count, paste0(output_prefix, "cell_count_prefilter.csv"), row.names = FALSE)

#### FILTERING ####

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 999999)
#seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200)

#### GENERATE POST-FILTER QC's ####
# VIOLIN PLOTS
post_filter_violin <- paste0(output_prefix, "violin_postfilter.pdf")

pdf(file = post_filter_violin,
    width = 8,
    height = 6)

plot(
  VlnPlot(object = seurat_obj,
          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
          ncol = 3,
          assay = "RNA")
)

dev.off()

# FEATURE SCATTER PLOTS
post_filter_scatter <- paste0(output_prefix, "scatter_postfilter.pdf")

pdf(file = post_filter_scatter,
    width = 8,
    height = 6)

plot1 <- FeatureScatter(object = seurat_obj,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.mt")

plot2 <- FeatureScatter(object = seurat_obj,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")

plot(plot1 + plot2) 

dev.off()

# DENSITY PLOTS - nFeature_RNA unscaled
post_filter_feature_dens <- paste0(output_prefix, "nFeature_unscaled_postfiltered_dens.pdf")

pdf(file = post_filter_feature_dens,
    width = 8,
    height = 6)

dens_plot_nfeature <- plotSingleCellDensity(input_obj = seurat_obj,
                                            metadata_variable = "nFeature_RNA",
                                            group.by = "orig.ident",
                                            scale_x_axis = FALSE) +
  geom_vline(xintercept = 200)

plot(dens_plot_nfeature)
dev.off()

# DENSITY PLOTS - nCount_RNA unscaled
post_filter_feature_dens <- paste0(output_prefix, "nCount_unscaled_postfiltered_dens.pdf")

pdf(file = post_filter_feature_dens,
    width = 8,
    height = 6)

dens_plot_ncount <- plotSingleCellDensity(input_obj = seurat_obj,
                                          metadata_variable = "nCount_RNA",
                                          group.by = "orig.ident",
                                          scale_x_axis = FALSE)  +
  geom_vline(xintercept = 500)

plot(dens_plot_ncount)
dev.off()

# DENSITY PLOTS - nFeature_RNA scaled
post_filter_feature_scaled_dens <- paste0(output_prefix, "nFeature_scaled_postfiltered_dens.pdf")

pdf(file = post_filter_feature_scaled_dens,
    width = 8,
    height = 6)

dens_plot_nfeature <- plotSingleCellDensity(input_obj = seurat_obj,
                                            metadata_variable = "nFeature_RNA",
                                            group.by = "orig.ident",
                                            scale_x_axis = TRUE) +
  geom_vline(xintercept = 200)

plot(dens_plot_nfeature)
dev.off()

# DENSITY PLOTS - nCount_RNA scaled
post_filter_feature_dens_scale <- paste0(output_prefix, "nCount_scaled_postfiltered_dens.pdf")

pdf(file = post_filter_feature_dens_scale,
    width = 8,
    height = 6)

dens_plot_ncount <- plotSingleCellDensity(input_obj = seurat_obj,
                                          metadata_variable = "nCount_RNA",
                                          group.by = "orig.ident",
                                          scale_x_axis = TRUE)  +
  geom_vline(xintercept = 500)

plot(dens_plot_ncount)
dev.off()

# Scatter Plot All 3 Variables

mt_postfilter <- FetchData(object = seurat_obj,
                          vars = c("orig.ident", 
                                   "nCount_RNA", 
                                   "nFeature_RNA", 
                                   "percent.mt")) %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point() + 
  scale_colour_gradient(low = "grey90", high = "black", limits = c(0,100)) +
  stat_smooth(method = lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident, scales = "free") +
  theme_classic()

mt_postfilter_name <- paste0(output_prefix, "mt_postfilter.pdf")
pdf(mt_postfilter_name,
    width = 8,
    height = 6)
plot(mt_postfilter)
dev.off()

# Save cell count prior to filtering
cell_count <- FetchData(object = seurat_obj, vars = "orig.ident") %>%
  tally(name = "Cells_after_filtering")

write.csv(cell_count, paste0(output_prefix, "cell_count_postfilter.csv"), row.names = FALSE)

#### NORMALIZATION ####
seurat_obj <- NormalizeData(seurat_obj, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)

#### FIND VARIABLE FEATURES ####
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst", 
                                   assay = "RNA",
                                   nfeatures = 2000)

#### SCALE DATA ####
# Scaling the data
all.genes <- rownames(seurat_obj)
# Default behavior is to only scale the 2000 most variable features
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

#### PCA #####
# WHY DIFFERENT?
seurat_obj <- RunPCA(seurat_obj, 
                     features = VariableFeatures(object = seurat_obj))
#seurat_obj <- RunPCA(seurat_obj, npcs = 50)

elbow_plot <- paste0(output_prefix, "/elbow.pdf")
pdf(elbow_plot,
    width = 8,
    height = 6)
plot(ElbowPlot(object = seurat_obj,
               ndims = 50,
               reduction = "pca"))
dev.off()

pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.03), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2


####

#### DIMENSIONS AND RESOLUTIONS ####

dim_output_dir <- paste0(output_prefix,"dims_res_explore/")
dir.create(dim_output_dir)

cluster_res <- seq(from = 0, to = 1, by = 0.1)
dims <- seq(from = 2, to = 15, by = 1)

seurat_obj_tmp <- seurat_obj

for (i in seq(cluster_res)){
  for (j in seq(dims)){
    
    seurat_obj_tmp <- FindNeighbors(seurat_obj_tmp,
                              dims = 1:dims[j])
    
    seurat_obj_tmp <- FindClusters(seurat_obj_tmp,
                                   resolution = cluster_res[i])

    seurat_obj_tmp <- RunUMAP(seurat_obj_tmp,
                              dims = 1:dims[j])
    
    file_name <- paste0(dim_output_dir,
                        "UMAP_dim_",
                        dims[j],
                        "_res_",
                        cluster_res[i],
                        ".png")
    
    png(file_name,
        width = 1200,
        height = 800)
    
    plot(DimPlot(seurat_obj_tmp,
                 label = TRUE,
                 reduction="umap") + 
          ggtitle(paste0("dim_res:", dims[j], "_", cluster_res[i])))
    
    dev.off()
  }
}

rm(seurat_obj_tmp)
gc()

#### CLUSTREE ####
dim_output_dir <- paste0(output_prefix,"dims_res_explore/")
dir.create(dim_output_dir)

tmp_seurat_clustree <- seurat_obj
clustree_res <- seq(from = 0, to = 1, by = 0.1)
dim <- 6
#dim <- 14

for (r in seq(clustree_res)){
  tmp_seurat_clustree <- RunUMAP(tmp_seurat_clustree, 
                                 reduction = "pca",
                                 dims = 1:dim)
  
  tmp_seurat_clustree <- FindNeighbors(tmp_seurat_clustree,
                                       reduction = "pca",
                                       dims = 1:dim)
  
  tmp_seurat_clustree <- FindClusters(tmp_seurat_clustree,
                                      resolution = clustree_res[r])
}

clustree_plot <- clustree(tmp_seurat_clustree)

pdf(paste0(dim_output_dir, "clustree.dim_", dim, ".pdf"),
    width = 8,
    height = 12)
plot(clustree_plot)
dev.off()

rm(tmp_seurat_clustree)
gc()

#### CLUSTER ####

#final_dim <- 5
final_res <- 0.7

final_dim <- 6
#final_res <- 0.4

seurat_obj <- RunUMAP(obj = seurat_obj,
                      reduction = "pca",
                      dims = 1:final_dim)

seurat_obj <- FindNeighbors(obj = seurat_obj,
                            reduction = "pca",
                            dims = 1:final_dim)

seurat_obj <- FindClusters(obj = seurat_obj,
                           resolution = final_res)

# Plot by various characteristics
#p1 <- DimPlot(seurat_obj,
#              group.by = "Group",
#              reduction = "umap")

p2 <- DimPlot(seurat_obj,
              group.by = "orig.ident",
              reduction = "umap")

p3 <- DimPlot(seurat_obj,
              reduction = "umap",
              label = TRUE)

png(paste0(output_prefix, "UMAPs_RNA.png"),
    width = 1600,
    height = 800)

#plot(p1 + p2 + p3)
plot(p2 + p3)
dev.off()

#### INITIAL CLUSTER QC ####
p1 <- DimPlot(seurat_obj,
              split.by = "orig.ident",
              reduction = "umap")

png(paste0(output_prefix, "UMAPS_RNA_by_sample.png"),
    width = 1800,
    height = 800)
plot(p1)
dev.off()

#### CELL MARKERS - FIND ALL MARKERS ####
seurat_obj.markers <- FindAllMarkers(seurat_obj,
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25)
# Gene markers across clusters

# Just looking at top 10
relevant_markers <- seurat_obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.csv(relevant_markers, paste0(output_prefix, "top_markers.csv"), row.names = FALSE)


#### MARKER FUNCTIONS ####

plot_markers_by_mtx <- function(seurat_obj_mtx,
                                marker_mtx,
                                feature,
                                out_dir_mtx = './marker_plots/', 
                                out_prefix_mtx = ''){
  # This will grab the correct columns from the matrix and try to plot the 
  # markers for the seurat object. This will also create help split up markers
  # in the case for the transcript so they are easily navigable.
  if (feature == 'gene'){
    plot_markers_by_list(seurat_obj_mtx,
                         unique(marker_mtx$gene),
                         paste0(out_dir_mtx, feature, '/'),
                         out_prefix_mtx)
    
  } else if (feature == 'transcript'){
    for (g in unique(marker_mtx$gene)){
      t <- marker_mtx[marker_mtx$gene == g,]$transcript
      
      # TODO: Find an R way to do this, need to find out how to check if markers exist
      tryCatch({
        plot_markers_by_list(seurat_obj_mtx,
                             unique(t),
                             paste0(out_dir_mtx, feature, '/'),
                             paste0(g, '.'))
      }, error = function(e){print(e)})
    }
  }
}

plot_markers_by_list <- function(seurat_obj_lst, 
                                 marker_lst,
                                 out_dir_lst = './marker_plots/',
                                 out_prefix_lst = ''){
  # Given a list of markers, this will generate a feature plot, violin plot,
  # and a dot plot for the seurat object that is passed in
  dir.create(file.path(out_dir_lst), recursive = TRUE)

  # Feature Plot
  png(paste0(out_dir_lst, out_prefix_lst, "markers_featureplot.png"),
      width = 1024,
      height = 768)
  
  print(
    FeaturePlot(object = seurat_obj_lst,
                features = marker_lst,
                reduction = "umap")
  )

  dev.off()
  
  # Violin Plot
  png(paste0(out_dir_lst, out_prefix_lst, "markers_violinplot.png"),
      width = 1024,
      height = 768)

  print(
    VlnPlot(object = seurat_obj_lst,
            split.by = "seurat_clusters",
            features = marker_lst,
            assay = "RNA")
  )

  dev.off()
  
  # Dot Plot
  png(paste0(out_dir_lst, out_prefix_lst, "markers_dotplot.png"),
      width = 1024,
      height = 768)
  
  print(
    DotPlot(object = seurat_obj_lst,
            features = marker_lst,
            group.by = "seurat_clusters",
            assay = "RNA",
            col.min= -1.5,
            col.max= 1.5,
            scale= TRUE,
            cluster.idents=FALSE,
            cols = "RdYlBu")  + coord_flip()
  )

  dev.off()
  
}

#### BLAZE MARKERS ####
markers <- read.table("input/references/blaze_markers.csv",
                      sep=",",
                      header = TRUE,
                      row.names = NULL)

plot_markers_by_mtx(seurat_obj,
                    markers,
                    feature_type,
                    paste0(output_prefix, "/markers/"),
                    '')
