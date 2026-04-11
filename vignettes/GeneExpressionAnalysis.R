## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(GeneExpressionAnalysis)
library(airway)

# Load the airway example data
data(airway)
expression_set <- loadfromSumE(airway)


## ----normalization------------------------------------------------------------
# View original distribution (using log_scale = TRUE to properly visualize raw counts)
plot_boxplot(expression_set, log_scale = TRUE)

# Remove outlier genes with a z-score threshold of 3
exp_clean <- remove_outliers(expression_set, threshold = 3)

# Filter out low expression genes, log transform, and normalize
exp_filtered <- filter_low_exp(exp_clean, min_count = 10)
exp_log <- log_transform(exp_filtered)
exp_norm <- quantile_norm(exp_log)

# View processed distribution (should be symmetrical and centered)
plot_boxplot(exp_norm, title = "Boxplot of Normalized Data")


## ----clustering---------------------------------------------------------------
# Plot the elbow to find optimal k
optimal_k(exp_norm, max_k = 10)

# Perform k-means clustering with 5 clusters
km <- kmeans_clust(exp_norm, k_clusters = 5)

# Perform hierarchical clustering and plot dendrogram
hc <- hierarchical_clust(exp_norm)
plot_den(hc$hclust)
plot_sample_den(exp_norm)


## ----pca----------------------------------------------------------------------
plot_PCA(exp_norm)


## ----enrichment, message=FALSE, warning=FALSE---------------------------------
# Extract the names of genes belonging to cluster 1
cluster1_genes <- names(km$cluster[km$cluster == 1])

# Perform GO enrichment analysis
go_res <- go_enrich(cluster1_genes)

# Plot the top enriched terms
plot_enrichment(go_res, top_n = 10)


## ----full-pipeline, eval=FALSE------------------------------------------------
# results <- run_full_pipeline(
#   expression_set = expression_set,
#   k_clusters = 5,
#   enrich_cluster = 1,
#   condition_col = "dex",
#   reference_level = "untrt"
# )
# 
# plot_heatmap(expression_set)
# plot_PCA(expression_set)

