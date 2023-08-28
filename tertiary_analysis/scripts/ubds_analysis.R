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

# The path to the input file
in_file <- "input/isoquant/q20.gene_counts.tsv"

# The name of the project to use in seurat
project_name <- "isoquant_q20_seurat"

# The prefix to output files
output_prefix <- "./output/isoquant/q20/gene/"

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


#### CELL MAKERS - BLAZE GENE MARKERS ####

markers <- read.table("input/references/blaze_markers.csv",
                      sep=",",
                      header = TRUE,
                      row.names = NULL)

gene_markers <- unique(markers$gene)

gene_markers_dir <- paste0(output_prefix, "/markers/")
dir.create(file.path(gene_marker_dir), recursive = TRUE)

# Feature Plot
png(paste0(gene_markers_dir, "blaze_gene_markers_featureplot.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = gene_markers,
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(gene_markers_dir, "blaze_gene_markers_violinplot.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = gene_markers,
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(gene_markers_dir, "blaze_gene_markers_dotplot.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = gene_markers,
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()

# TODO: Refine the transcript markers stuff to go off the marker list

#### CELL MAKERS - BLAZE TRANSCRIPT MARKERS (ENSG00000007372) ####
markers.blaze <-c("ENST00000606377.7",
                  "ENST00000639386.2",
                  "ENST00000419022.6",
                  "ENST00000643871.1",
                  "ENST00000638914.3",
                  "ENST00000640610.1",
                  "ENST00000379132.8",
                  "ENST00000379109.7",
                  "ENST00000379129.7",
                  "ENST00000639916.1",
                  "ENST00000640368.2",
                  "ENST00000638685.1",
                  "ENST00000379107.7",
                  "ENST00000639409.1",
                  "ENST00000640975.1",
                  "ENST00000638963.1",
                  "ENST00000638965.1",
                  "ENST00000241001.13",
                  "ENST00000474783.2",
                  "ENST00000638629.1",
                  "ENST00000639548.1",
                  "ENST00000640125.1",
                  "ENST00000481563.6",
                  "ENST00000638696.1",
                  "ENST00000638755.1",
                  "ENST00000470027.7",
                  "ENST00000640287.1",
                  "ENST00000533333.5",
                  "ENST00000640613.1",
                  "ENST00000423822.7",
                  "ENST00000639061.1",
                  "ENST00000640766.1",
                  "ENST00000640963.1",
                  "ENST00000471303.6",
                  "ENST00000379111.7",
                  "ENST00000639950.1",
                  "ENST00000464174.6",
                  "ENST00000494377.7",
                  "ENST00000638250.1",
                  "ENST00000379123.10",
                  "ENST00000640038.1",
                  "ENST00000524853.6",
                  "ENST00000640872.1",
                  "ENST00000639394.1",
                  "ENST00000640735.1",
                  "ENST00000639109.1",
                  "ENST00000638853.1",
                  "ENST00000640460.1",
                  "ENST00000638913.1",
                  "ENST00000638802.1",
                  "ENST00000638346.1",
                  "ENST00000530373.6",
                  "ENST00000379115.9",
                  "ENST00000639006.1",
                  "ENST00000639054.1",
                  "ENST00000640335.1",
                  "ENST00000638762.1",
                  "ENST00000438681.6",
                  "ENST00000638878.1",
                  "ENST00000639943.1",
                  "ENST00000640431.1",
                  "ENST00000532916.5",
                  "ENST00000639034.2",
                  "ENST00000531910.6",
                  "ENST00000640172.1",
                  "ENST00000640684.1",
                  "ENST00000639079.1",
                  "ENST00000640242.1",
                  "ENST00000455099.6",
                  "ENST00000530714.6",
                  "ENST00000534353.5",
                  "ENST00000533156.2",
                  "ENST00000639920.1",
                  "ENST00000534390.2",
                  "ENST00000527769.5",
                  "ENST00000532175.5",
                  "ENST00000640251.1",
                  "ENST00000638278.1",
                  "ENST00000640617.1",
                  "ENST00000639203.1",
                  "ENST00000640819.1")

# Feature Plot
png(paste0(output_prefix, "blaze_transcript_markers_featureplot.ENSG00000007372.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = unique(markers.blaze),
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(output_prefix, "blaze_transcript_markers_violinplot.ENSG00000007372.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = unique(markers.blaze),
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(output_prefix, "blaze_transcript_markers_dotplot.ENSG00000007372.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = unique(markers.blaze),
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()

#### CELL MAKERS - BLAZE TRANSCRIPT MARKERS (ENSG00000026025) ####
markers.blaze <-c("ENST00000544301.7",
                  "ENST00000478746.1",
                  "ENST00000497849.1",
                  "ENST00000224237.9",
                  "ENST00000487938.5",
                  "ENST00000485947.1",
                  "ENST00000469543.5",
                  "ENST00000421459.2",
                  "ENST00000637053.1",
                  "ENST00000495528.1")

# Feature Plot
png(paste0(output_prefix, "blaze_transcript_markers_featureplot.ENSG00000026025.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = unique(markers.blaze),
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(output_prefix, "blaze_transcript_markers_violinplot.ENSG00000026025.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = unique(markers.blaze),
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(output_prefix, "blaze_transcript_markers_dotplot.ENSG00000026025.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = unique(markers.blaze),
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()
#### CELL MAKERS - BLAZE TRANSCRIPT MARKERS (ENSG00000104888) ####
markers.blaze <-c("ENST00000221485.8",
                  "ENST00000600601.5",
                  "ENST00000600672.5",
                  "ENST00000596689.1",
                  "ENST00000598018.1")

# Feature Plot
png(paste0(output_prefix, "blaze_transcript_markers_featureplot.ENSG00000104888.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = unique(markers.blaze),
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(output_prefix, "blaze_transcript_markers_violinplot.ENSG00000104888.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = unique(markers.blaze),
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(output_prefix, "blaze_transcript_markers_dotplot.ENSG00000104888.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = unique(markers.blaze),
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()
#### CELL MAKERS - BLAZE TRANSCRIPT MARKERS (ENSG00000124785) ####
markers.blaze <-c("ENST00000244766.7",
                  "ENST00000622188.4",
                  "ENST00000616243.1",
                  "ENST00000495850.1")

# Feature Plot
png(paste0(output_prefix, "blaze_transcript_markers_featureplot.ENSG00000124785.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = unique(markers.blaze),
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(output_prefix, "blaze_transcript_markers_violinplot.ENSG00000124785.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = unique(markers.blaze),
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(output_prefix, "blaze_transcript_markers_dotplot.ENSG00000124785.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = unique(markers.blaze),
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()
#### CELL MAKERS - BLAZE TRANSCRIPT MARKERS (ENSG00000131095) ####

# These markers are the transcripts of the gene in the header + all overlapping
# transcripts. These were included because when comparing only the transcripts,
# the markers were not found so the marker list was expanded to include all
# overlapping markers
#markers.blaze <-c("ENST00000417826.3",
#                  "ENST00000588735.3",
#                  "ENST00000639277.1",
#                  "ENST00000638304.1",
#                  "ENST00000441312.2",
#                  "ENST00000639243.1",
#                  "ENST00000638488.1",
#                  "ENST00000639369.1",
#                  "ENST00000592065.2",
#                  "ENST00000640859.1",
#                  "ENST00000638400.1",
#                  "ENST00000639042.1",
#                  "ENST00000640545.1",
#                  "ENST00000253408.11",
#                  "ENST00000638921.1",
#                  "ENST00000586125.2",
#                  "ENST00000638618.1",
#                  "ENST00000589701.2",
#                  "ENST00000592706.5",
#                  "ENST00000585543.6",
#                  "ENST00000591880.2",
#                  "ENST00000588640.5",
#                  "ENST00000640552.1",
#                  "ENST00000435360.8",
#                  "ENST00000591327.2",
#                  "ENST00000638281.1",
#                  "ENST00000639921.1",
#                  "ENST00000592320.6",
#                  "ENST00000586127.6",
#                  "ENST00000586793.6",
#                  "ENST00000587997.6",
#                  "ENST00000376990.8",
#                  "ENST00000591719.5",
#                  "ENST00000590922.1",
#                  "ENST00000588957.5",
#                  "ENST00000588316.1",
#                  "ENST00000585728.5",
#                  "ENST00000588037.1",
#                  "ENST00000593179.1",
#                  "ENST00000331733.5",
#                  "ENST00000577339.5",
#                  "ENST00000410006.6",
#                  "ENST00000357776.6",
#                  "ENST00000417826.3",
#                  "ENST00000410027.5",
#                  "ENST00000331733.5",
#                  "ENST00000426333.7",
#                  "ENST00000592576.5",
#                  "ENST00000591382.5",
#                  "ENST00000591382.5",
#                  "ENST00000588374.1",
#                  "ENST00000593072.5",
#                  "ENST00000589825.5",
#                  "ENST00000592408.5",
#                  "ENST00000592701.2",
#                  "ENST00000590105.1",
#                  "ENST00000587309.5",
#                  "ENST00000593135.6",
#                  "ENST00000590129.1")

# This marker list contains ONLY the transcripts that are apart of the gene
markers.blaze <-c("ENST00000417826.3",
                  "ENST00000588735.3",
                  "ENST00000639277.1",
                  "ENST00000638304.1",
                  "ENST00000441312.2",
                  "ENST00000639243.1",
                  "ENST00000638488.1",
                  "ENST00000639369.1",
                  "ENST00000592065.2",
                  "ENST00000640859.1",
                  "ENST00000638400.1",
                  "ENST00000639042.1",
                  "ENST00000640545.1",
                  "ENST00000253408.11",
                  "ENST00000638921.1",
                  "ENST00000586125.2",
                  "ENST00000638618.1",
                  "ENST00000589701.2",
                  "ENST00000592706.5",
                  "ENST00000585543.6",
                  "ENST00000591880.2",
                  "ENST00000588640.5",
                  "ENST00000640552.1",
                  "ENST00000435360.8",
                  "ENST00000591327.2",
                  "ENST00000638281.1",
                  "ENST00000639921.1",
                  "ENST00000592320.6",
                  "ENST00000586127.6",
                  "ENST00000586793.6",
                  "ENST00000587997.6",
                  "ENST00000376990.8",
                  "ENST00000591719.5",
                  "ENST00000590922.1",
                  "ENST00000588957.5",
                  "ENST00000588316.1",
                  "ENST00000585728.5",
                  "ENST00000588037.1",
                  "ENST00000593179.1")

# Feature Plot
png(paste0(output_prefix, "blaze_transcript_markers_featureplot.ENSG00000131095.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = unique(markers.blaze),
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(output_prefix, "blaze_transcript_markers_violinplot.ENSG00000131095.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = unique(markers.blaze),
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(output_prefix, "blaze_transcript_markers_dotplot.ENSG00000131095.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = unique(markers.blaze),
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()
#### CELL MAKERS - BLAZE TRANSCRIPT MARKERS (ENSG00000176165) ####

# These markers are the transcripts of the gene in the header + all overlapping
# transcripts. These were included because when comparing only the transcripts,
# the markers were not found so the marker list was expanded to include all
# overlapping markers
#markers.blaze <-c("ENST00000706482.1",
#                  "ENST00000313071.7",
#                  "ENST00000658593.1",
#                  "ENST00000551395.5",
#                  "ENST00000546560.5",
#                  "ENST00000549487.1",
#                  "ENST00000706482.1",
#                  "ENST00000313071.7",
#                  "ENST00000675861.1",
#                  "ENST00000675861.1",
#                  "ENST00000622740.3",
#                  "ENST00000399387.9",
#                  "ENST00000668749.1",
#                  "ENST00000689292.1",
#                  "ENST00000548213.3",
#                  "ENST00000552957.6",
#                  "ENST00000653638.1",
#                  "ENST00000671672.1")

# This marker list contains ONLY the transcripts that are apart of the gene
markers.blaze <-c("ENST00000706482.1",
                  "ENST00000313071.7")

# Feature Plot
png(paste0(output_prefix, "blaze_transcript_markers_featureplot.ENSG00000176165.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = unique(markers.blaze),
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(output_prefix, "blaze_transcript_markers_violinplot.ENSG00000176165.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = unique(markers.blaze),
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(output_prefix, "blaze_transcript_markers_dotplot.ENSG00000176165.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = unique(markers.blaze),
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()
#### CELL MAKERS - BLAZE TRANSCRIPT MARKERS (ENSG00000162374) ####
markers.blaze <-c("ENST00000651693.1",
                  "ENST00000463650.2",
                  "ENST00000448907.7",
                  "ENST00000652252.1",
                  "ENST00000371827.5",
                  "ENST00000651347.1",
                  "ENST00000357083.8",
                  "ENST00000650764.1",
                  "ENST00000494555.2",
                  "ENST00000371824.7",
                  "ENST00000371823.8",
                  "ENST00000652693.1",
                  "ENST00000371819.1",
                  "ENST00000651258.1",
                  "ENST00000652353.1",
                  "ENST00000371821.6",
                  "ENST00000652274.1",
                  "ENST00000492299.2",
                  "ENST00000474675.1")

# Feature Plot
png(paste0(output_prefix, "blaze_transcript_markers_featureplot.ENSG00000162374.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = unique(markers.blaze),
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(output_prefix, "blaze_transcript_markers_violinplot.ENSG00000162374.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = unique(markers.blaze),
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(output_prefix, "blaze_transcript_markers_dotplot.ENSG00000162374.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = unique(markers.blaze),
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()
#### CELL MAKERS - BLAZE TRANSCRIPT MARKERS (ENSG00000171786) ####
markers.blaze <-c("ENST00000302101.6")

# Feature Plot
png(paste0(output_prefix, "blaze_transcript_markers_featureplot.ENSG00000171786.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = unique(markers.blaze),
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(output_prefix, "blaze_transcript_markers_violinplot.ENSG00000171786.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = unique(markers.blaze),
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(output_prefix, "blaze_transcript_markers_dotplot.ENSG00000171786.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = unique(markers.blaze),
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()
#### CELL MAKERS - BLAZE TRANSCRIPT MARKERS (ENSG00000148123) ####
markers.blaze <-c("ENST00000374874.8",
                  "ENST00000456287.5",
                  "ENST00000494890.5",
                  "ENST00000395056.2",
                  "ENST00000463206.1")

# Feature Plot
png(paste0(output_prefix, "blaze_transcript_markers_featureplot.ENSG00000148123.png"), width = 1024, height = 768)

FeaturePlot(object = seurat_obj,
            features = unique(markers.blaze),
            reduction = "umap")

dev.off()

# Violin Plot
png(paste0(output_prefix, "blaze_transcript_markers_violinplot.ENSG00000148123.png"), width = 1024, height = 768)

VlnPlot(object = seurat_obj,
        split.by = "seurat_clusters",
        features = unique(markers.blaze),
        assay = "RNA")

dev.off()

# Dot Plot
png(paste0(output_prefix, "blaze_transcript_markers_dotplot.ENSG00000148123.png"), width = 1024, height = 768)

DotPlot(object = seurat_obj,
        features = unique(markers.blaze),
        group.by = "seurat_clusters",
        assay = "RNA",
        col.min= -1.5,
        col.max= 1.5,
        scale= TRUE,
        cluster.idents=FALSE,
        cols = "RdYlBu")  + coord_flip()

dev.off()