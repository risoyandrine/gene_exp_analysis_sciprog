library(GeneExpressionAnalysis)
library(Biobase)

# we create a small mock ExpressionSet for testing
set.seed(42)
mock_data <- matrix(runif(100, 10, 1000), nrow = 10)
rownames(mock_data) <- paste0("gene", 1:10)
colnames(mock_data) <- paste0("sample", 1:10)

# we inject a row of pure zeros to test filtering
mock_data[2, ] <- 0

mock_pheno <- data.frame(
  condition = rep(c("treated", "untreated"), each = 5),
  row.names = colnames(mock_data)
)
mock_eset <- ExpressionSet(assayData = mock_data, phenoData = AnnotatedDataFrame(mock_pheno))

test_that("filter_low_exp drops low expression genes", {
  res <- filter_low_exp(mock_eset, min_count = 5)
  expect_s4_class(res, "ExpressionSet")
  # gene 2 was all zeros, so it should be filtered out
  expect_false("gene2" %in% rownames(exprs(res)))
  # should have 9 genes left
  expect_equal(nrow(exprs(res)), 9)
})

test_that("log_transform correctly scales values", {
  res <- log_transform(mock_eset, pseudo_count = 1)
  expect_s4_class(res, "ExpressionSet")

  # ensure the values were actually log scaled (max should be much smaller than 1000)
  expect_true(max(exprs(res)) < 50)

  # spot check math: log2(0 + 1) = 0
  expect_equal(exprs(res)["gene2", "sample1"], 0)
})

test_that("clustering module outputs correct structures", {
  # for clustering to work well, we need data
  exp_norm <- log_transform(mock_eset)

  # check k-means
  km <- kmeans_clust(exp_norm, k_clusters = 2)
  expect_s3_class(km, "kmeans")
  expect_length(km$cluster, 10)

  # check hierarchical
  hc <- hierarchical_clust(exp_norm)
  expect_type(hc, "list")
  expect_s3_class(hc$hclust, "hclust")
  expect_s3_class(hc$dist, "dist")
})

test_that("remove_outliers executes without destroying structure", {
  res <- remove_outliers(mock_eset, threshold = 1) # strict threshold
  expect_s4_class(res, "ExpressionSet")
  # since it removes genes, the number of features should be less than or equal to original
  expect_true(nrow(exprs(res)) <= nrow(exprs(mock_eset)))
})

test_that("quantile_norm equalizes sample distributions", {
  res <- quantile_norm(log_transform(filter_low_exp(mock_eset)))
  medians <- apply(exprs(res), 2, median)
  expect_true(diff(range(medians)) < 1e-10)
})

test_that("expression_summary returns invisible NULL and prints correctly", {
  output <- capture.output(res <- expression_summary(mock_eset), type = "output")
  expect_null(res)
  # check that the output contains expected labels
  expect_true(any(grepl("Number of genes", output)))
  expect_true(any(grepl("Number of samples", output)))
  expect_true(any(grepl("Mean expression", output)))
})
