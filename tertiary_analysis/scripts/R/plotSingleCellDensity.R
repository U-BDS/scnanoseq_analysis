#' plotSingleCellDensity
#' This generic helper function creates a density plot of the distribution of cells (nFeature) or molecules/UMIs (nCount)
#' from an input object which may be a data.frame with single-cell/nuclei metadata or a Seurat object.
#' @param input_obj Input object (data.frame of single-cell metadata or Seurat obj.)
#' @param metadata_variable Metadata column name linked to cell numbers and UMIs - typically `nFeature_assay_name` or `nCount_assay_name`
#' @param group.by Metadata column name linked to the grouping variable to group the plot by
#' @param scale_x_axis Transforms x-axis to the log10 scale
#' @param geom_density_level Set transparency level for density plot - lower values generate more transparent density curves
#'
#' @return A density plot
#' @export
#'
#' @examples
#' \dontrun{
#' plotSingleCellDensity(input_obj, metadata_variable = "nFeature_RNA", group.by = "Stim")
#' }
plotSingleCellDensity <- function(input_obj, ...) {
  UseMethod("plotSingleCellDensity")
}

plotSingleCellDensity.Seurat <- function(input_obj,
                                         metadata_variable,
                                         group.by,
                                         scale_x_axis = FALSE,
                                         geom_density_level = 0.2) {
  metadata <- dplyr::select(
    input_obj@meta.data,
    {{ metadata_variable }},
    {{ group.by }}
  )

  meta_density <- ggplot2::ggplot(
    metadata,
    aes(
      x = .data[[metadata_variable]],
      color = .data[[group.by]],
      fill = .data[[group.by]]
    )
  ) +
    geom_density(alpha = geom_density_level) +
    theme_classic()

  if (scale_x_axis == TRUE) {
    return(meta_density + scale_x_log10())
  } else {
    return(meta_density)
  }
}

plotSingleCellDensity.data.frame <- function(input_obj,
                                             metadata_variable,
                                             group.by,
                                             scale_x_axis = FALSE,
                                             geom_density_level = 0.2) {
  metadata <- dplyr::select(
    input_obj,
    {{ metadata_variable }},
    {{ group.by }}
  )

  meta_density <- ggplot2::ggplot(
    metadata,
    aes(
      x = .data[[metadata_variable]],
      color = .data[[group.by]],
      fill = .data[[group.by]]
    )
  ) +
    geom_density(alpha = geom_density_level) +
    theme_classic()

  if (scale_x_axis == TRUE) {
    return(meta_density + scale_x_log10())
  } else {
    return(meta_density)
  }
}
