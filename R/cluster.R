#' Performs a Kmeans clustering on the expression data
#'
#' @param expression_set An expression set object
#' @param k_clusters Number of clusters, default is set to 5
#' @param seed Seed for reproducibility, default is set to 42
#' @return A kmeans object
#' @examples
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' km_result <- kmeans_clust(expression_set, 5)
#' km_clust <- km_result$cluster
#' @export

kmeans_clust <- function(expression_set, k_clusters = 5, seed = 42) {
  exp <- Biobase::exprs(expression_set)
  # we start with some basic checks to validate input
  if (!is.numeric(exp)) stop("Expression data must be a numeric matrix")
  if (nrow(exp) == 0) stop("No genes remaining after filtering, we cannot perform clustering")
  if (k_clusters > nrow(exp)) stop("Number of clusters is greater than the number of genes")

  set.seed(seed)

  result <- stats::kmeans(exp,
    centers = k_clusters,
    iter.max = 100,
    nstart = 10
  )
  return(result)
}


#' Performs a hierarchical clustering on the expression data
#'
#' @param expression_set An expression set object
#' @param method The hierarchical clustering method to be used, set to complete as default
#' @return A hierarchical clustering object
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' hc_tree <- hierarchical_clust(expression_set, "complete")
#' hc_clust <- cutree(hc_tree, k = 5)
#' }
#' @export

hierarchical_clust <- function(expression_set, method = "complete") {
  exp <- Biobase::exprs(expression_set)

  distance_matrix <- dist(exp)
  hierarchicalclust <- hclust(distance_matrix, method = method)

  return(list(hclust = hierarchicalclust, dist = distance_matrix))
}

#' Plot the Elbow Method to find the optimal number of clusters
#'
#' @param expression_set An expression set object
#' @param max_k The maximum number of clusters to test, default is 10
#' @param seed Seed for reproducibility
#' @return A ggplot object showing the elbow plot
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' optimal_k(expression_set)
#' }
#' @export
optimal_k <- function(expression_set, max_k = 10, seed = 42) {
  exp <- Biobase::exprs(expression_set)
  set.seed(seed)
  # we calculate Total Within-Cluster Sum of Squares (WSS)
  wss <- numeric(max_k)
  for (k in 1:max_k) {
    km_result <- stats::kmeans(exp, centers = k, iter.max = 100, nstart = 10)
    wss[k] <- km_result$tot.withinss
  }
  plot_data <- data.frame(k = 1:max_k, wss = wss)
  elbow_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = k, y = wss)) +
    ggplot2::geom_line(color = "slateblue3", linewidth = 1) +
    ggplot2::geom_point(color = "coral4", size = 3) +
    ggplot2::scale_x_continuous(breaks = 1:max_k) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Elbow Plot for Choosing Number of Clusters",
      x = "Number of Clusters",
      y = "Total Within-Cluster Sum of Squares"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(elbow_plot)
}
