#' This function will load gene expression data from a CSV file and read in ExpressionSet
#'
#' @param file Represents a path to a CSV file with genes as rows and samples as columns
#' @return An ExpressionSet object
#' @examples
#' \donttest{
#' expression_set <- loadfromCSV("path/to/data.csv")
#' }
#' @export
loadfromCSV <- function(file) {
  if (!file.exists(file)) {
    stop("File was not found: ", file)
  }
  expression_matrix <- as.matrix(read.csv(file, row.names = 1)) # read the CSV file
  expression_set <- Biobase::ExpressionSet(assayData = expression_matrix)

  return(expression_set)
}

#' This function can read a summarized experiment object (like we are using as an example, the airway dataset) and convert to expression set
#'
#' @param SumE A SummarizedExperiment object
#' @return An ExpressionSet object
#' @examples
#' data(example_airway)
#' expression_set <- example_airway
#' @export

loadfromSumE <- function(SumE) {
  if (!is(SumE, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object")
  }
  expression_matrix <- SummarizedExperiment::assay(SumE)
  pheno_data <- as.data.frame(SummarizedExperiment::colData(SumE)) # we extract sample information
  pheno <- Biobase::AnnotatedDataFrame(data = pheno_data)
  expression_set <- Biobase::ExpressionSet(assayData = expression_matrix, phenoData = pheno)

  return(expression_set)
}
