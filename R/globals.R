#' @importFrom stats dist hclust prcomp kmeans cutree var
#' @importFrom utils read.csv head
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#' @importFrom tidyr pivot_longer
NULL

utils::globalVariables(c("pc_1", "pc_2", "gene", "expression", "sample"))
