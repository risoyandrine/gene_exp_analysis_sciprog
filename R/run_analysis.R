#' The function runs the full pipeline
#'
#' @param expression_set An expression set object
#' @param k_clusters Number of clusters, default is set to 5
#' @param seed Seed for reproducibility, default is set to 42
#' @param method The hierarchical clustering method to be used, default is set to "complete"
#' @param top_n Number of top enriched terms to plot, default is set to 10
#' @param enrich_cluster The cluster to be used for the enrichment analysis, default is set to 1
#' @param condition_col The column in the phenoData to be used for the enrichment analysis, default is set to "dex"
#' @param reference_level The reference level for the enrichment analysis, default is set to "untrt"
#' @param OrgDb The organism database to be used for the enrichment analysis, default is set to "org.Hs.eg.db"
#' @param organism The organism to be used for the enrichment analysis, default is set to "hsa"
#' @param keyType Type of gene ID, default is "ENSEMBL"
#' @param already_log Logical, whether the data is already log-transformed. Default is FALSE.
#' @param min_count minimum count threshold for filtering, default is set to 10
#' @param min_samples minimum number of samples for filtering, default is set to 2
#' @param top_n_genes Number of genes to show in heatmap, default is 1000
#' @param top_genes_pca Number of genes to use for PCA, default is 500
#' @param threshold Standard deviation threshold for outlier detection, default is 3
#' @return A list containing:
#' \itemize{
#'   \item expression_set: The processed ExpressionSet
#'   \item kmeans: The kmeans clustering result
#'   \item hclust: The hierarchical clustering result (list with hclust and dist)
#'   \item go: GO enrichment results
#'   \item kegg: KEGG enrichment results
#'   \item gsea: GSEA enrichment results
#' }
#' @examples
#' if (requireNamespace("airway", quietly = TRUE)) {
#'   library(airway)
#'   data(airway)
#'   expression_set <- loadfromSumE(airway)
#'   run_full_pipeline(expression_set)
#' }
#' @export


run_full_pipeline <- function(expression_set, k_clusters = 5,
                              seed = 42, method = "complete", top_n = 10,
                              enrich_cluster = 1, condition_col = "dex",
                              reference_level = "untrt",
                              OrgDb = "org.Hs.eg.db",
                              organism = "hsa", keyType = "ENSEMBL",
                              already_log = FALSE, min_count = 10,
                              min_samples = 2, top_n_genes = 1000,
                              top_genes_pca = 500, threshold = 3) {
  # safety check
  if (enrich_cluster > k_clusters) {
    stop("enrich_cluster cannot be greater than k_clusters")
  }

  # step 1: preprocessing of the data
  expression_set <- remove_outliers(expression_set, threshold = threshold)
  expression_set <- filter_low_exp(
    expression_set,
    min_count = min_count,
    min_samples = min_samples,
    already_log = already_log
  )

  if (!already_log) {
    expression_set <- log_transform(expression_set)
  }
  expression_set <- quantile_norm(expression_set)

  # step 2: clustering
  optimal_k(expression_set, max_k = 10)
  km <- kmeans_clust(expression_set, k_clusters, seed)
  hc <- hierarchical_clust(expression_set, method)

  # step 3: enrichment analysis
  go_result <- go_enrich(names(km$cluster[km$cluster == enrich_cluster]),
    OrgDb = OrgDb, keyType = keyType
  )
  kegg_result <- kegg_enrich(names(km$cluster[km$cluster == enrich_cluster]),
    OrgDb = OrgDb,
    organism = organism,
    keyType = keyType
  )
  gse_result <- go_gse(expression_set, condition_col, reference_level,
    OrgDb = OrgDb, keyType = keyType
  )
  # step 4: plotting
  print(plot_distr(expression_set))
  print(plot_boxplot(expression_set, title = "Boxplot of Normalized Data"))
  plot_heatmap(expression_set, top_n_genes = top_n_genes)
  print(plot_PCA(expression_set, top_genes_pca = top_genes_pca))

  gene_var <- apply(Biobase::exprs(expression_set), 1, var)
  exp_top <- expression_set[order(gene_var, decreasing = TRUE)[seq_len(50)], ]
  hc_vis <- hierarchical_clust(exp_top, method)
  plot_den(hc_vis$hclust)

  plot_sample_den(expression_set)
  print(plot_enrichment(go_result, top_n = top_n))
  print(plot_enrichment(kegg_result, top_n = top_n))
  print(plot_enrichment(gse_result, top_n = top_n))
  cat("\nExploratory Analysis Summary:\n")
  print(expression_summary(expression_set))
  cat("\nClustering Analysis Summary:\n")
  print(table(km$cluster))

  cat("\n Enrichment Analysis\n")

  print_enrich_summary <- function(result_obj, name) {
    df <- as.data.frame(result_obj)
    cat(sprintf("Significant %s terms: %d\n", name, nrow(df)))
    if (nrow(df) > 0) {
      cat("Top hits:", paste(head(df$Description, 3), collapse = " | "), "\n")
    }
  }

  print_enrich_summary(go_result, "GO")
  print_enrich_summary(kegg_result, "KEGG")
  print_enrich_summary(gse_result, "GSEA")
  return(list(
    expression_set = expression_set,
    kmeans = km,
    hclust = hc,
    go = go_result,
    kegg = kegg_result,
    gsea = gse_result
  ))
}
