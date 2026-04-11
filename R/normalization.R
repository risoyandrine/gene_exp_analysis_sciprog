#' Normalize the data, we want to filter low expression genes from the ExpressionSet
#'
#' @param expression_set Expression Set object
#' @param min_count minimum count threshold, default is set to 10
#' @param min_samples minimum number of samples that a gene have to be expressed in, default is set to 2
#' @param already_log Logical, whether the data is already log-transformed. Default is FALSE.
#' @return a filtered expression set object without the lowly expressed genes
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' exp_set_filtered <- filter_low_exp(expression_set, already_log = FALSE)
#' }
#' @export


filter_low_exp <- function(expression_set, min_count = 10, min_samples = 2, already_log = FALSE) {
  exp <- Biobase::exprs(expression_set)
  # we start with some basic checks
  if (!is.numeric(Biobase::exprs(expression_set))) stop("Expression data must be numeric")
  if (nrow(Biobase::exprs(expression_set)) == 0) stop("No genes remaining after filtering")

  # if data is already log-transformed, we adjust the min_count threshold proportionally
  if (already_log) {
    message("Data is log-transformed. Adjusting filtering threshold.")
    min_count <- log2(min_count + 1)
  }

  genes_to_keep <- rowSums(exp >= min_count) >= min_samples
  message(sum(!genes_to_keep), " genes were removed, ", sum(genes_to_keep), " genes were kept") # how many genes were kept
  return(expression_set[genes_to_keep, ])
}

#' Next step is to logtransform the data
#'
#' @param expression_set Expression Set object
#' @param pseudo_count to avoid log(0)
#' @return a log transformed expression set object
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' exp_set_filtered <- filter_low_exp(expression_set)
#' exp_set_log <- log_transform(exp_set_filtered)
#' }
#' @export


log_transform <- function(expression_set, pseudo_count = 1) {
  Biobase::exprs(expression_set) <- log2(Biobase::exprs(expression_set) + pseudo_count)

  return(expression_set)
}

#' Lastly, for the preprocessing of the data, we perform a quantile normalization
#'
#' @param expression_set Expression Set object
#' @return a quantile normalized expression set object
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' exp_set_filtered <- filter_low_exp(expression_set)
#' exp_set_log <- log_transform(exp_set_filtered)
#' exp_set_norm <- quantile_norm(exp_set_log)
#' }
#' @export

quantile_norm <- function(expression_set) {
  Biobase::exprs(expression_set) <- limma::normalizeBetweenArrays(Biobase::exprs(expression_set), method = "quantile")
  return(expression_set)
}
