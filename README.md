# GeneExpressionAnalysis

`GeneExpressionAnalysis` is an R package designed for gene expression analysis. It provides a pipeline for preprocessing, clustering, and biological interpretation of RNA-seq data within the Bioconductor ecosystem.

The package simplifies workflows by combining outlier detection, normalization, multiple clustering methods (K-means and Hierarchical), and downstream enrichment analysis (GO, KEGG, and GSEA) into a suite of tools.

## Features

*   **Exploratory Data Analysis (EDA)**: Summary statistics and distribution plots (density and boxplots).
*   **Outlier Detection**: Identification and removal of outlier genes and samples using Z-score thresholds.
*   **Data Normalization**: Low-expression filtering, log-transformation, and quantile normalization for cross-sample comparability.
*   **Clustering**: 
    *   **K-means**: Includes elbow plot analysis (`optimal_k`) to determine the best cluster number.
    *   **Hierarchical Clustering**: Supports various linkage methods for both genes and samples.
*   **Biological Enrichment**: 
    *   **GO Enrichment**: Over-representation analysis (ORA) for Biological Processes.
    *   **KEGG Enrichment**: Pathway analysis with automated gene ID conversion.
    *   **GSEA**: Gene Set Enrichment Analysis to identify subtle, coordinated changes.
*   **Visualization**: Includes PCA plots, heatmaps, and dendrograms.

## Installation

Install the latest version of `GeneExpressionAnalysis` from GitHub:

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("andrinerisoy/GeneExpressionAnalysis")
```

## Quick Start

The fastest way to analyze your data is using the `run_full_pipeline` function, which executes the entire workflow and generates a summary of plots and results.

```R
library(GeneExpressionAnalysis)
library(airway)

# Load example data
data(airway)
expression_set <- loadfromSumE(airway)

# Execute the full pipeline
results <- run_full_pipeline(
  expression_set = expression_set, 
  k_clusters = 5, 
  enrich_cluster = 1,
  condition_col = "dex", 
  reference_level = "untrt"
)

# Access all results in a single list
normalized_data <- results$expression_set
kmeans_res      <- results$kmeans
hclust_res      <- results$hclust
go_results      <- results$go
kegg_results    <- results$kegg
gsea_results    <- results$gsea
```

## Detailed Workflow

For maximum control, you can run each step of the analysis manually:

### 1. Exploratory Analysis & Outlier Detection
```R
# View data summary
expression_summary(expression_set)

# Visualize distributions
plot_boxplot(expression_set, log_scale = TRUE)
plot_distr(expression_set)

# Identify outliers
gene_outliers <- gene_outlier_detection(expression_set, threshold = 3)
sample_outliers <- sample_outlier_detection(expression_set, threshold = 3)

# Filter outlier, only genes are removed, samples are flagged
exp_clean <- remove_outliers(expression_set, threshold = 3)
```

### 2. Normalization
```R
# Filter low-expression, log-transform, and quantile normalize
exp_filtered <- filter_low_exp(exp_clean, min_count = 10)
exp_log <- log_transform(exp_filtered)
exp_norm <- quantile_norm(exp_log)
```

### 3. Clustering
```R
# K-means with Elbow Plot
optimal_k(exp_norm, max_k = 10)
km <- kmeans_clust(exp_norm, k_clusters = 5)

# Hierarchical Clustering
hc <- hierarchical_clust(exp_norm, method = "complete")
```

### 4. Visualization
```R
# Principal Component Analysis
plot_PCA(exp_norm)

# Heatmap of top 50 most variable genes
plot_heatmap(exp_norm, top_n_genes = 50)

# Dendrograms
plot_den(hc$hclust) # Gene dendrogram
plot_sample_den(exp_norm) # Sample dendrogram
```

### 5. Enrichment Analysis
```R
# Extract genes from a cluster of interest
cluster_genes <- names(km$cluster[km$cluster == 1])

# Run specific enrichment methods
go_res <- go_enrich(cluster_genes)
kegg_res <- kegg_enrich(cluster_genes)
gsea_res <- go_gse(exp_norm, condition_col = "dex", reference_level = "untrt")

# Plot results
plot_enrichment(go_res)
```

---
**Author**: Andrine Risøy  
**License**: GPL-3  
**Project**: Scientific Programming in R 
**Date**: April 2026