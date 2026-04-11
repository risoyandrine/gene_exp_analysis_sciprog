#' This function plots the distribution of expression values per sample
#'
#' @param expression_set An expression set object
#' @param title Plot title for the distribution plot
#' @return a ggplot2 plot of the expression distribution
#' @examples
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' expression_set <- filter_low_exp(expression_set)
#' expression_set <- log_transform(expression_set)
#' plot_distr(expression_set)
#' @export


plot_distr <- function(expression_set, title = "Distribution of Expression Data") {
  exp <- Biobase::exprs(expression_set)

  # we want to plot with ggplot2 and therefore need to transform the data to correct format
  exp_df <- as.data.frame(exp)
  exp_df$gene <- rownames(exp_df)
  exp_long <- tidyr::pivot_longer(exp_df, -gene, names_to = "sample", values_to = "expression")

  ggplot2::ggplot(exp_long, ggplot2::aes(x = expression, fill = sample)) +
    ggplot2::geom_density(alpha = 0.3) +
    ggplot2::labs(title = title, x = "Expression", y = "Density") +
    ggplot2::facet_wrap(~sample) +
    ggplot2::theme_minimal()
}


#' This function will detect potential outliers from the data
#'
#' @param expression_set An expression set object
#' @param threshold Tolerated standard deviations from the mean
#' @return Outlier genes
#' @examples
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' expression_set <- filter_low_exp(expression_set)
#' expression_set <- log_transform(expression_set)
#' gene_outlier_detection(expression_set, threshold = 3)
#' @export


gene_outlier_detection <- function(expression_set, threshold = 3) {
  exp <- Biobase::exprs(expression_set)

  # calculate mean and sd for each gene
  avg_gene_exp <- rowMeans(exp, na.rm = TRUE)

  # calculate z-scores across genes
  z_scores <- scale(avg_gene_exp)
  z_scores <- as.numeric(z_scores)

  # identify outliers above set threshold
  outliers <- rownames(exp)[abs(z_scores) > threshold]

  message(length(outliers), " outlier genes detected with Z-score bigger than ", threshold)
  return(outliers)
}

#' This function will detect high variance samples
#'
#' @param expression_set An expression set object
#' @param threshold Tolerated standard deviations from the mean
#' @return Outlier samples
#' @examples
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' expression_set <- filter_low_exp(expression_set)
#' expression_set <- log_transform(expression_set)
#' sample_outlier_detection(expression_set, threshold = 3)
#' @export


sample_outlier_detection <- function(expression_set, threshold = 3) {
  exp <- Biobase::exprs(expression_set)

  # calculate mean and sd for each sample
  avg_sample_exp <- colMeans(exp, na.rm = TRUE)

  # calculate z-scores across genes
  z_scores <- scale(avg_sample_exp)
  z_scores <- as.numeric(z_scores)

  # identify outliers above set threshold
  outliers <- colnames(exp)[abs(z_scores) > threshold]

  message(length(outliers), " outlier samples detected with Z-score bigger than ", threshold)
  return(outliers)
}

#' Summary of the expression data
#'
#' @param expression_set An expression set object
#' @return Invisible NULL, prints the summary of the expression data to console
#' @examples
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' expression_summary(expression_set)
#' @export


expression_summary <- function(expression_set) {
  exp <- Biobase::exprs(expression_set)

  cat("Summary of the Expression data:\n")
  cat("Number of genes:   ", nrow(exp), "\n")
  cat("Number of samples: ", ncol(exp), "\n")
  cat("Min expression:    ", round(min(exp, na.rm = TRUE), 2), "\n")
  cat("Max expression:    ", round(max(exp, na.rm = TRUE), 2), "\n")
  cat("Mean expression:   ", round(mean(exp, na.rm = TRUE), 2), "\n")
  cat("Sample names:      ", colnames(exp), "\n")

  invisible(NULL) # for a cleaner output
}

#' Remove outlier genes from the expression set
#'
#' @param expression_set An expression set object
#' @param threshold Tolerated standard deviations from the mean
#' @return A filtered expression set object with outlier genes removed
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' expression_set <- filter_low_exp(expression_set)
#' exp_clean <- remove_outliers(expression_set, threshold = 3)
#' }
#' @export
remove_outliers <- function(expression_set, threshold = 3) {
  gene_outliers <- gene_outlier_detection(expression_set, threshold = threshold)
  sample_outliers <- sample_outlier_detection(expression_set, threshold = threshold)

  if (length(sample_outliers) > 0) {
    # we keep only samples that are not in the outlier list
    warning("Outlier samples detected, but not removed. Consider removing outliers manually: ", paste(sample_outliers, collapse = ", "))
  }

  if (length(gene_outliers) > 0) {
    # we keep only genes that are not in the outlier list
    genes_to_keep <- !(rownames(Biobase::exprs(expression_set)) %in% gene_outliers)
    expression_set <- expression_set[genes_to_keep, ]
    message("Removed ", length(gene_outliers), " outlier genes.")
  } else {
    message("No outliers to remove.")
  }

  return(expression_set)
}

#' A function to plot the expression distribution using boxplots per sample
#'
#' @param expression_set An expression set object
#' @param title Plot title for the boxplot
#' @param log_scale Logical, whether to apply a log2 transformation before plotting. Default is FALSE.
#' @return a ggplot2 boxplot of the expression distribution
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' expression_set <- filter_low_exp(expression_set)
#' plot_boxplot(expression_set, log_scale = TRUE)
#' }
#' @export
plot_boxplot <- function(expression_set, title = "Boxplot of Expression Data", log_scale = FALSE) {
  exp <- Biobase::exprs(expression_set)

  # we want to plot with ggplot2 and therefore need to transform the data to correct format
  exp_df <- as.data.frame(exp)
  exp_df$gene <- rownames(exp_df)
  exp_long <- tidyr::pivot_longer(exp_df, -gene, names_to = "sample", values_to = "expression")

  if (log_scale) {
    exp_long$expression <- log2(exp_long$expression + 1)
    y_label <- "Log2 Expression"
  } else {
    y_label <- "Expression"
  }

  ggplot2::ggplot(exp_long, ggplot2::aes(x = sample, y = expression, fill = sample)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.2) +
    ggplot2::labs(title = title, x = "Sample", y = y_label) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}
