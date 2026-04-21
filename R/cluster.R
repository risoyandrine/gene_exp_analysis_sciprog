#' Performs a Kmeans clustering on the expression data
#'
#' @param expression_set An expression set object
#' @param k_clusters Number of clusters, default is set to 5
#' @return A kmeans object
#' @examples
#' data(example_airway)
#' expression_set <- example_airway
#' km_result <- kmeans_clust(expression_set, 5)
#' km_clust <- km_result$cluster
#' @export

kmeans_clust <- function(expression_set, k_clusters = 5) {
  exp <- Biobase::exprs(expression_set)
  # check that the input is valid
  if (!is.numeric(exp)) stop("The expression data has to be a numeric matrix")
  if (nrow(exp) == 0) stop("No genes are remaining after filtering")
  if (k_clusters > nrow(exp)) stop("The number of clusters is larger than the number of genes")

  result <- stats::kmeans(exp,
    centers = k_clusters,
    iter.max = 100,
    nstart = 25
  )
  return(result)
}


#' Performs a hierarchical clustering on the expression data
#'
#' @param expression_set An expression set object
#' @param method The hierarchical clustering method to be used, set to complete as default
#' @return A hierarchical clustering object
#' @examples
#' data(example_airway)
#' expression_set <- example_airway
#' expression_set <- filter_low_exp(expression_set)
#' expression_set <- log_transform(expression_set)
#' expression_set <- quantile_norm(expression_set)
#' gene_var <- apply(Biobase::exprs(expression_set), 1, var)
#' expression_set <- expression_set[order(gene_var, decreasing = TRUE)[1:50], ]
#' hc_tree <- hierarchical_clust(expression_set, "complete")
#' hc_clust <- cutree(hc_tree$hclust, k = 5)
#' @export

hierarchical_clust <- function(expression_set, method = "complete") {
  if (!is(expression_set, "ExpressionSet")) stop("The expression data need to be an ExpressionSet object")
  exp <- Biobase::exprs(expression_set)

  distance_matrix <- dist(exp)
  hierarchicalclust <- hclust(distance_matrix, method = method)

  return(list(hclust = hierarchicalclust, dist = distance_matrix))
}

#' Plot the Elbow Method to find the optimal number of clusters
#'
#' @param expression_set An expression set object
#' @param max_k The maximum number of clusters to test, default is 10
#' @return A ggplot object showing the elbow plot
#' @examples
#' data(example_airway)
#' expression_set <- example_airway
#' optimal_k(expression_set)
#' @export
optimal_k <- function(expression_set, max_k = 10) {
  exp <- Biobase::exprs(expression_set)
  wss <- numeric(max_k)
  for (k in seq_len(max_k)) {
    km_result <- stats::kmeans(exp, centers = k, iter.max = 100, nstart = 25)
    wss[k] <- km_result$tot.withinss
  }
  plot_data <- data.frame(k = seq_len(max_k), wss = wss)
  elbow_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = k, y = wss)) +
    ggplot2::geom_line(color = "slateblue3", linewidth = 1) +
    ggplot2::geom_point(color = "coral4", size = 3) +
    ggplot2::scale_x_continuous(breaks = seq_len(max_k)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Elbow Plot for Choosing Number of Clusters",
      x = "Number of Clusters",
      y = "Total Within-Cluster Sum of Squares"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(elbow_plot)
}
