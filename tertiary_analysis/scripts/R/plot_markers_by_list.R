#' plot_markers_by_list
#' Plots a series of figures for marker visualization: FeaturePlot, violin plots and dot plot
#' @param seurat_obj seurat object
#' @param marker_lst vector of features to plot
#' @param out_dir_lst output directory
#' @param out_prefix_lst output prefix
#'
#' @return
#' @export

plot_markers_by_list <- function(seurat_obj, 
                                 marker_lst,
                                 out_dir_lst = './marker_plots/',
                                 out_prefix_lst = ''){
  # Given a list of markers, this will generate a feature plot, violin plot,
  # and a dot plot for the seurat object that is passed in
  dir.create(file.path(out_dir_lst), recursive = TRUE, showWarnings = FALSE)
  
  # Feature Plot
  png(paste0(out_dir_lst, out_prefix_lst, "markers_featureplot.png"),
      width = 2000,
      height = 1000)
  
  plot(
    FeaturePlot(object = seurat_obj,
                features = marker_lst,
                reduction = "umap")
  )
  
  dev.off()
  
  # Violin Plot
  png(paste0(out_dir_lst, out_prefix_lst, "markers_violinplot.png"),
      width = 2000,
      height = 1000)
  
  plot(
    VlnPlot(object = seurat_obj,
            split.by = "seurat_clusters",
            features = marker_lst,
            assay = "RNA")
  )
  
  dev.off()
  
  # Dot Plot
  png(paste0(out_dir_lst, out_prefix_lst, "markers_dotplot.png"),
      width = 2000,
      height = 1000)
  
  plot(
    DotPlot(object = seurat_obj,
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