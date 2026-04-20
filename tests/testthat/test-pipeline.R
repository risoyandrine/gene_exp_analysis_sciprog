library(GeneExpressionAnalysis)
library(Biobase)

# mock ExpressionSet for testing
mock_data <- matrix(runif(100, 10, 1000), nrow = 10)
rownames(mock_data) <- paste0("gene", 1:10)
colnames(mock_data) <- paste0("sample", 1:10)

# gene2 is all zeros so we can verify filtering
mock_data[2, ] <- 0

mock_pheno <- data.frame(
  condition = rep(c("treated", "untreated"), each = 5),
  row.names = colnames(mock_data)
)
mock_eset <- ExpressionSet(assayData = mock_data, phenoData = AnnotatedDataFrame(mock_pheno))

# loading of the data

test_that("loadfromCSV loads a CSV into an ExpressionSet", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(mock_data, file = tmp)

  res <- loadfromCSV(tmp)
  expect_s4_class(res, "ExpressionSet")
  expect_equal(nrow(exprs(res)), nrow(mock_data))
  expect_equal(ncol(exprs(res)), ncol(mock_data))

  unlink(tmp)
})

test_that("loadfromCSV errors on missing file", {
  expect_error(loadfromCSV("nonexistent.csv"), "File was not found")
})

# Preprocessing of the data including filtering, log transform, normalization, outliers

test_that("filter_low_exp drops low expression genes", {
  res <- filter_low_exp(mock_eset, min_count = 5)
  expect_s4_class(res, "ExpressionSet")
  expect_false("gene2" %in% rownames(exprs(res)))
  expect_equal(nrow(exprs(res)), 9)
})

test_that("log_transform correctly scales values", {
  res <- log_transform(mock_eset, pseudo_count = 1)
  expect_s4_class(res, "ExpressionSet")
  expect_true(max(exprs(res)) < 50)
  # log2(0 + 1) = 0
  expect_equal(exprs(res)["gene2", "sample1"], 0)
})

test_that("quantile_norm equalizes sample distributions", {
  res <- quantile_norm(log_transform(filter_low_exp(mock_eset)))
  medians <- apply(exprs(res), 2, median)
  expect_true(diff(range(medians)) < 1e-10)
})

test_that("gene_outlier_detection returns gene names", {
  res <- gene_outlier_detection(mock_eset, threshold = 0.5)
  expect_type(res, "character")
  expect_true(all(res %in% rownames(exprs(mock_eset))))
})

test_that("sample_outlier_detection returns sample names", {
  res <- sample_outlier_detection(mock_eset, threshold = 0.5)
  expect_type(res, "character")
  expect_true(all(res %in% colnames(exprs(mock_eset))))
})

test_that("remove_outliers preserves ExpressionSet structure", {
  res <- remove_outliers(mock_eset, threshold = 1)
  expect_s4_class(res, "ExpressionSet")
  expect_true(nrow(exprs(res)) <= nrow(exprs(mock_eset)))
})

# testing of the exploratory analysis on the data

test_that("expression_summary prints and returns invisible NULL", {
  output <- capture.output(res <- expression_summary(mock_eset), type = "output")
  expect_null(res)
  expect_true(any(grepl("Number of genes", output)))
  expect_true(any(grepl("Number of samples", output)))
  expect_true(any(grepl("Mean expression", output)))
})

# testing of the clustering functions on the data

test_that("kmeans_clust returns a valid kmeans object", {
  exp_norm <- log_transform(mock_eset)
  km <- kmeans_clust(exp_norm, k_clusters = 2)
  expect_s3_class(km, "kmeans")
  expect_length(km$cluster, 10)
})

test_that("hierarchical_clust returns hclust and dist", {
  exp_norm <- log_transform(mock_eset)
  hc <- hierarchical_clust(exp_norm)
  expect_type(hc, "list")
  expect_s3_class(hc$hclust, "hclust")
  expect_s3_class(hc$dist, "dist")
})

test_that("optimal_k returns an elbow plot", {
  exp_norm <- log_transform(mock_eset)
  p <- optimal_k(exp_norm, max_k = 3)
  expect_s3_class(p, "ggplot")
})

# testing of the visualization functions on the data

test_that("plot_distr returns a ggplot", {
  p <- plot_distr(log_transform(mock_eset))
  expect_s3_class(p, "ggplot")
})

test_that("plot_boxplot returns a ggplot", {
  p <- plot_boxplot(mock_eset)
  expect_s3_class(p, "ggplot")
})

test_that("plot_boxplot works with log_scale", {
  p <- plot_boxplot(mock_eset, log_scale = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap runs without error", {
  exp_norm <- quantile_norm(log_transform(filter_low_exp(mock_eset, min_count = 5)))
  expect_no_error(plot_heatmap(exp_norm, top_n_genes = 5))
})

test_that("plot_PCA returns a ggplot", {
  exp_norm <- quantile_norm(log_transform(filter_low_exp(mock_eset, min_count = 5)))
  p <- plot_PCA(exp_norm, top_genes_pca = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_den returns an hclust object invisibly", {
  hc <- hierarchical_clust(log_transform(mock_eset))
  res <- plot_den(hc$hclust)
  expect_s3_class(res, "hclust")
})

test_that("plot_sample_den returns an hclust object invisibly", {
  res <- plot_sample_den(log_transform(mock_eset))
  expect_s3_class(res, "hclust")
})

# testing of the full pipeline which is the most efficient way to run the analysis

test_that("run_full_pipeline executes end-to-end", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("airway")

  library(airway)
  data(airway)
  eset <- loadfromSumE(airway)

  result <- run_full_pipeline(
    eset,
    k_clusters      = 3,
    condition_col   = "dex",
    reference_level = "untrt",
    top_n_genes     = 500,
    top_genes_pca   = 200
  )

  expect_type(result, "list")
  expect_named(result, c("expression_set", "kmeans", "hclust", "go", "kegg", "gsea"),
    ignore.order = TRUE
  )
  expect_s4_class(result$expression_set, "ExpressionSet")
  expect_s3_class(result$kmeans, "kmeans")
  expect_equal(length(unique(result$kmeans$cluster)), 3)
  expect_s3_class(result$hclust$hclust, "hclust")
})

test_that("run_full_pipeline rejects invalid enrich_cluster", {
  expect_error(
    run_full_pipeline(mock_eset, k_clusters = 2, enrich_cluster = 5),
    "enrich_cluster cannot be greater than k_clusters"
  )
})
