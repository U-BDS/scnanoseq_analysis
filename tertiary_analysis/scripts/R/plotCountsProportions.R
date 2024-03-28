#' plotCountsProportions
#' This functions fetches counts data from a Seurat object and plots the proportion of counts relative to the selected metadata.
#' @param input_obj Input object
#' @param metadata_variable Metadata column name linked to plot proportions relative to (e.g: cell name or anatomical regions)
#' @param group.by Metadata column name linked to the grouping variable to group the plot by (orig.ident or experimental group)
#'
#' @return A bar chart plot
#' @export
#'
#' @examples
#' \dontrun{
#' plotCountsProportions(input_obj, metadata_variable = "seurat_clusters", group.by = "Stim")
#' }
plotCountsProportions <- function(input_obj,
                                  metadata_variable,
                                  group.by) {
  
  metadata <- FetchData(input_obj, vars = c( group.by, metadata_variable )) %>% 
    group_by( .data[[group.by]], .data[[metadata_variable]] ) %>% 
    summarize(total = n())  
  
  # compute total cell/spot numbers
  metadata %>%
    group_by( .data[[metadata_variable]] ) %>%
    summarise(total_by_meta = sum(total)) -> total_by_meta
  
  # add totals in metadata data.frame
  metadata <- left_join(metadata, total_by_meta, metadata_variable)
  
  # compute proportions
  metadata$proportions <- (metadata$total/metadata$total_by_meta)
  
  proportion_plot <- ggplot(metadata, aes(x = .data[[metadata_variable]], y = proportions, fill = .data[[group.by]], label = proportions)) +
    geom_col(color = "black") +
    geom_text(label = round(metadata$proportions, digits = 2), position = position_stack(vjust = .5)) +
    theme_bw() 
  
  return(proportion_plot)
  
}
