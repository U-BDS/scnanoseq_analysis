#' plotCounts
#' This generic function creates a bar chart of the number of cells or spots grouped by a variable
#' from the input object which may be a data.frame with single-cell/nuclei/spatial metadata or a Seurat object.
#' @param input_obj Input object (data.frame of single-cell metadata or Seurat obj.)
#' @param metadata_variable Metadata column name linked to cell names, anatomical region or cluster names
#' @param group.by Metadata column name linked to the grouping variable to group the plot by
#'
#' @return A bar chart plot
#' @export
#'
#' @examples
#' \dontrun{
#' plotCounts(input_obj, metadata_variable = "seurat_clusters", group.by = "Stim")
#' }
plotCounts <- function(input_obj, ...) {
  UseMethod("plotCounts")
}

plotCounts.Seurat <- function(input_obj,
                              metadata_variable,
                              group.by) {
  metadata <- dplyr::select(
    input_obj@meta.data,
    {{ metadata_variable }},
    {{ group.by }}
  )
  
  meta_histogram <- ggplot2::ggplot(
    metadata,
    aes(x = .data[[metadata_variable]], fill = .data[[group.by]])
  ) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
  
  return(meta_histogram)
}

plotCounts.data.frame <- function(input_obj,
                                  metadata_variable,
                                  group.by) {
  metadata <- dplyr::select(
    input_obj,
    {{ metadata_variable }},
    {{ group.by }}
  )
  
  meta_histogram <- ggplot2::ggplot(
    metadata,
    aes(x = .data[[metadata_variable]], fill = .data[[group.by]])
  ) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
  
  return(meta_histogram)
}