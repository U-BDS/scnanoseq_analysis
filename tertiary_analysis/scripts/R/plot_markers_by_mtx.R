#' plot_markers_by_mtx
#' For gene or transcript feature, this functions calls plot_markers_by_list
#' @param seurat_obj seurat object
#' @param marker_mtx data.frame of gene and transcript markers. A columns named  "gene_id", "transcript_id" are expected to be present
#' @param feature feature type: `gene` or `transcript` (or `trans`)
#' @param out_dir_mtx output directory
#' @param out_prefix_mtx output prefix
#'
#' @return
#' @export

plot_markers_by_mtx <- function(seurat_obj,
                                marker_mtx,
                                feature,
                                out_dir_mtx = './marker_plots/', 
                                out_prefix_mtx = ''){
  # This will grab the correct columns from the matrix and try to plot the 
  # markers for the seurat object. This will also create help split up markers
  # in the case for the transcript so they are easily navigable.
  if (feature == 'gene'){
    
    # filter marker to only IDs present in obj:
    status <- unique(marker_mtx$gene_id) %in% rownames(seurat_obj@assays$RNA@counts)
    markers_filtered <- unique(marker_mtx$gene_id)[status]
    
    if (length(markers_filtered > 0)) {
      
      plot_markers_by_list(seurat_obj,
                           markers_filtered,
                           paste0(out_dir_mtx, feature, '/'),
                           out_prefix_mtx)
      
    } else {
      
      message("No genes from input markers were found in the current object")
      
    }
    
  } else if (feature == 'transcript' | feature == 'trans') {
    
     for (g in unique(marker_mtx$gene_id)){
       transcript_by_gene <- marker_mtx[marker_mtx$gene_id == g,]$transcript_id
       
       # filter marker to only IDs present in obj:
       status <- unique(transcript_by_gene) %in% rownames(seurat_obj@assays$RNA@counts)
       markers_filtered <- unique(transcript_by_gene)[status]
       
       if (length(markers_filtered > 0)) {
         
         plot_markers_by_list(seurat_obj,
                              markers_filtered,
                              paste0(out_dir_mtx, feature, '/'),
                              paste0(g, '.'))
         
       } else {
         
         message(paste0("Transcripts from gene ", g, " were not found in the current object"))
         
       }
    }
  }
}
