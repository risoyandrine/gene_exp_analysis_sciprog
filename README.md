# GeneExpressionAnalysis

This is my R project for Scientific Programming. The project involves the package`GeneExpressionAnalysis` for performing gene expression analysis. The package aim to provide a pipeline for preprocessing, clustering, and biological interpretation of RNA-seq data.

Overall, it simplifies workflows by combining outlier detection, normalization, clustering methods (K-means and Hierarchical), and downstream enrichment analysis (GO, KEGG, and GSEA) to a single package.

## Features

*   **Exploratory Data Analysis (EDA)**: Performs summary statistics and distribution plots to get some insights of the data.
*   **Outlier Detection**: Identification and removal of outlier genes and samples using Z-score thresholds, specifically outlier samples are flagged and outlier genes are removed.
*   **Data Normalization**: Low-expression filtering, log-transformation, and quantile normalization to make samples comparable for downstream analysis.
*   **Clustering**: 
    *   **K-means**: Includes elbow plot analysis (`optimal_k`) to determine the ideal number of clusters.
    *   **Hierarchical Clustering**: Supports different linkage methods for both genes and samples.
*   **Biological Enrichment**: 
    *   **GO Enrichment**: Over-representation analysis (ORA) for Biological Processes.
    *   **KEGG Enrichment**: Pathway analysis with gene ID conversion.
    *   **GSEA**: Gene Set Enrichment Analysis to identify subtle, coordinated changes in gene expression.
*   **Visualization**: Includes PCA plots, heatmaps, and dendrograms.

## Installation

Below is the code to install the latest version of `GeneExpressionAnalysis` from GitHub:

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("risoyandrine/gene_exp_analysis_sciprog")```

## Quick Start

To make the pipeline more efficient, the `run_full_pipeline` function was created. This function executes the entire workflow and generates a summary of plots and results.

```R
library(GeneExpressionAnalysis)
library(airway)

# Load the example data
data(airway)
expression_set <- loadfromSumE(airway)

# An example for running the full pipeline with some default parameters
results <- run_full_pipeline(
  expression_set = expression_set, 
  k_clusters = 5, 
  enrich_cluster = 1,
  condition_col = "dex", 
  reference_level = "untrt"
)

# The results are stored in a single list
normalized_data <- results$expression_set
kmeans_res      <- results$kmeans
hclust_res      <- results$hclust
go_results      <- results$go
kegg_results    <- results$kegg
gsea_results    <- results$gsea
```

## Step-by-step workflow
If you want to run the pipeline step-by-step, the workflow is as follows:

1. Exploratory Data Analysis
2. Outlier Detection
3. Data Normalization
4. Clustering
5. Biological Enrichment
6. Visualization


---
**Author**: Andrine Risøy  
**Date**: April 2026