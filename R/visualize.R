#' Plot heatmap using pheatmap
#'
#' @param expression_set An expression set object
#' @param show_rownames logical, whether to show row names
#' @param cluster_rows logical, whether to cluster rows
#' @param cluster_cols logical, whether to cluster columns
#' @param scale scale, whether to scale the data
#' @param top_n_genes Integer, the number of top variable genes to plot if the dataset is large. Default is 1000.
#' @return a heatmap of the gene expression
#' @examples
#' \donttest{
#' data(example_airway, package = "GeneExpressionAnalysis")
#' expression_set <- example_airway
#' exp_set_filtered <- filter_low_exp(expression_set)
#' exp_set_log <- log_transform(exp_set_filtered)
#' exp_set_norm <- quantile_norm(exp_set_log)
#' plot_heatmap(exp_set_norm)
#' }
#' @export


plot_heatmap <- function(expression_set,
                         show_rownames = FALSE,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         scale = "row",
                         top_n_genes = 1000) {
  exp <- Biobase::exprs(expression_set)

  gene_vars <- apply(exp, 1, var, na.rm = TRUE)
  exp <- exp[gene_vars > 0 & !is.na(gene_vars), ]

  # prevent memory exhaustion on large datasets by visualizing only the top variable genes
  if (nrow(exp) > top_n_genes) {
    message(sprintf(
      "Large dataset detected (> %d genes),
                    subsetting to top %d most variable genes for heatmap visualization.",
      top_n_genes, top_n_genes
    ))
    gene_vars_current <- apply(exp, 1, var, na.rm = TRUE)
    exp <- exp[order(gene_vars_current, decreasing = TRUE)[seq_len(top_n_genes)], ]
  }

  # custom blue to white to red color palette
  my_colors <- colorRampPalette(c("navy", "white", "red"))(50)

  pheatmap::pheatmap(exp,
    show_rownames = show_rownames,
    cluster_rows = cluster_rows, # cluster genes
    cluster_cols = cluster_cols, # cluster samples
    scale = scale, # scale the data, row sets the mean to 0 and sd to 1 for each gene
    color = my_colors,
    main = "Heatmap of gene expression"
  )
}

#' Plot PCA plot
#'
#' @param expression_set An expression set object
#' @param top_genes_pca The number of top variable genes to use for PCA. Default is 500.
#' @return a PCA plot of the gene expression
#' @examples
#' data(example_airway, package = "GeneExpressionAnalysis")
#' expression_set <- example_airway
#' exp_set_filtered <- filter_low_exp(expression_set)
#' exp_set_log <- log_transform(exp_set_filtered)
#' exp_set_norm <- quantile_norm(exp_set_log)
#' plot_PCA(exp_set_norm)
#' @export


plot_PCA <- function(expression_set, top_genes_pca = 500) {
  exp <- Biobase::exprs(expression_set)

  gene_vars <- apply(exp, 1, var)
  exp <- exp[gene_vars > 0, ]

  if (nrow(exp) > top_genes_pca) {
    exp <- exp[order(gene_vars[gene_vars > 0],
      decreasing = TRUE
    )[seq_len(top_genes_pca)], ]
  }

  pca <- prcomp(t(exp), scale. = TRUE)

  explained_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  df_pca <- data.frame(
    pc_1 = pca$x[, 1],
    pc_2 = pca$x[, 2],
    sample = rownames(pca$x)
  )


  ggplot2::ggplot(
    df_pca,
    ggplot2::aes(x = pc_1, y = pc_2, label = sample, color = sample)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -1, show.legend = FALSE) +
    ggplot2::labs(
      x = paste0("PC1: ", explained_var[1], "% Var Exp"),
      y = paste0("PC2: ", explained_var[2], "% Var Exp")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA plot")
}


#' Plot dendrogram
#'
#' @param hc An hclust object
#' @return plots the dendrogram
#' @examples
#' data(example_airway, package = "GeneExpressionAnalysis")
#' expression_set <- example_airway
#' exp_set_filtered <- filter_low_exp(expression_set)
#' exp_set_log <- log_transform(exp_set_filtered)
#' exp_set_norm <- quantile_norm(exp_set_log)
#' hc <- hierarchical_clust(exp_set_norm)
#' plot_den(hc$hclust)
#' @export

plot_den <- function(hc) {
  plot(hc, main = "Gene clustering dendrogram")
  invisible(hc)
}

#' Plot dendrogram of samples for better interpretation of the data
#'
#' @param expression_set An expression set object
#' @param method The hierarchical clustering method to be used, default is "complete"
#' @return plots the dendrogram
#' @examples
#' data(example_airway, package = "GeneExpressionAnalysis")
#' expression_set <- example_airway
#' exp_set_filtered <- filter_low_exp(expression_set)
#' exp_set_log <- log_transform(exp_set_filtered)
#' exp_set_norm <- quantile_norm(exp_set_log)
#' plot_sample_den(exp_set_norm)
#' @export

plot_sample_den <- function(expression_set, method = "complete") {
  exp <- Biobase::exprs(expression_set)
  hc <- hclust(dist(t(exp)), method = method) # the t() is for transposing the matrix to cluster samples instead of genes
  plot(hc, main = "Sample clustering dendrogram", sub = "", xlab = "")
  invisible(hc)
}
