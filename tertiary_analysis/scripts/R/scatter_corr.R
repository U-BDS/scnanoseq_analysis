#' scatter_corr
#' Creates a scatter plot with linear regression and Pearson correlation value between 2 features
#' @param data.frame data.frame contains features of interest
#' @param feature feature to be plotted on x-axis
#' @param feature2 feature to be plotted on y-axis
#' @param feature_type feature type ("gene" or "transcript" level) for visualization purposes
#' @param axis_size size of font in axis
#' @param title_size size of font in title
#'
#' @return
#' @export
#'
#' @examples
scatter_corr <- function(data.frame, feature, feature2,
                         feature_type = c("gene", "transcript"),
                         axis_size = 20, title_size = 24) {
  
  # set visuals based on feature_type
  if (feature_type == "gene") {
    color <- "#0072B2"
    title <- "Gene-level"
  } else {
    color <- "#009E73"
    title <- "Transcript-level"
  }
  
  #plot
  scatter_plot <- ggplot(data.frame,
                         aes(x=.data[[feature]],
                             y=.data[[feature2]])) +
    geom_point(size = 3, alpha = 0.7, stroke = 1, colour = color) +
    stat_smooth(method = lm, colour = "#383434") +
    theme(axis.text.x = element_text(size=axis_size),
          axis.text.y = element_text(size=axis_size),
          axis.title =  element_text(size=title_size),
          axis.line = element_line(),
          plot.title = element_text(size = title_size, face = "bold"),
          panel.background = element_blank(),
          legend.key = element_blank()) +
    ggtitle(
      paste0(title, "\nPearson Correlation: ",
             round(cor(data.frame[,feature],
                       data.frame[,feature2],
                       method = "pearson"), digits = 2))
    )
  
  return(scatter_plot)
}