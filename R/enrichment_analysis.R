#' Performs a GO enrichment analysis
#'
#' @param gene_list A list of genes
#' @param OrgDb OrgDb database to use for annotation, default is "org.Hs.eg.db"
#' @param keyType Type of gene ID, default is "ENSEMBL"
#' @return A GO enrichment object
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' km <- kmeans_clust(expression_set, 5)
#' gene_list <- names(km$cluster[km$cluster == 1])
#' go_enrich(gene_list)
#' }
#' @export

go_enrich <- function(gene_list, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL") {
  enrichment_go <- clusterProfiler::enrichGO(gene_list, OrgDb = OrgDb, ont = "BP", keyType = keyType)
  return(enrichment_go)
}

#' Performs a KEGG enrichment analysis
#'
#' @param gene_list A list of genes
#' @param OrgDb OrgDb database to use for annotation, default is "org.Hs.eg.db"
#' @param organism KEGG organism code, default is "hsa"
#' @param keyType Type of gene ID, default is "ENSEMBL"
#' @return A KEGG enrichment object
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' km <- kmeans_clust(expression_set, 5)
#' gene_list <- names(km$cluster[km$cluster == 1])
#' kegg_enrich(gene_list)
#' }
#' @export

kegg_enrich <- function(gene_list, OrgDb = "org.Hs.eg.db", organism = "hsa", keyType = "ENSEMBL") {
  entrez_ids <- clusterProfiler::bitr(gene_list, fromType = keyType, toType = "ENTREZID", OrgDb = OrgDb) # convert to entrez ids, which are KEGG compatible format
  if (nrow(entrez_ids) == 0) stop("No genes could be converted to Entrez IDs")
  enrichment_kegg <- clusterProfiler::enrichKEGG(gene = entrez_ids$ENTREZID, organism = organism)
  return(enrichment_kegg)
}


#' Performs Gene Set Enrichment Analysis
#'
#' @param expression_set An expression set object
#' @param condition_col The column in the phenoData will be used for the differential expression
#' @param reference_level The reference level for the condition
#' @param OrgDb OrgDb database to use for annotation, default is "org.Hs.eg.db"
#' @param keyType Type of gene ID, default is "ENSEMBL"
#' @return A GSEA enrichment object
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' exp_norm <- quantile_norm(log_transform(filter_low_exp(expression_set)))
#' go_gse(exp_norm, condition_col = "dex", reference_level = "untrt")
#' }
#' @export

go_gse <- function(expression_set, condition_col, reference_level, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL") {
  exp <- Biobase::exprs(expression_set)
  groups <- Biobase::pData(expression_set)[[condition_col]]

  group_means <- tapply(seq_len(ncol(exp)), groups, function(i) rowMeans(exp[, i, drop = FALSE]))

  treat_level <- setdiff(unique(groups), reference_level)
  log2fc <- group_means[[treat_level]] - group_means[[reference_level]]

  ranked <- sort(log2fc, decreasing = TRUE)

  result <- clusterProfiler::gseGO(ranked, ont = "BP", OrgDb = OrgDb, keyType = keyType)
  return(result)
}

#' Plot the enrichment results
#'
#' @param enrichment_result An enrichment result object
#' @param top_n Number of top enriched terms to plot
#' @return A plot of the enrichment results
#' @examples
#' \dontrun{
#' library(airway)
#' data(airway)
#' expression_set <- loadfromSumE(airway)
#' km <- kmeans_clust(expression_set, 5)
#' gene_list <- names(km$cluster[km$cluster == 1])
#' go_result <- go_enrich(gene_list)
#' plot_enrichment(go_result)
#' }
#' @export

plot_enrichment <- function(enrichment_result, top_n = 10) {
  # if no significant pathways were found, return an error message
  if (nrow(as.data.frame(enrichment_result)) == 0) {
    message("No significant enrichment pathways found to plot.")
    return(invisible(NULL))
  }

  enrichplot::dotplot(enrichment_result, showCategory = top_n) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Enrichment analysis")
}
