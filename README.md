# Illumina-Microarray-Analysis-Pipeline

This repository provides an R script designed for analyzing Illumina gene expression data retrieved from the Gene Expression Omnibus (GEO). The script includes robust data preprocessing, statistical analysis, and visualization tools.

## Features
- **Data Retrieval**: Automatically fetches datasets from GEO using GEOquery.
- **Preprocessing**: Handles data normalization, transformation, and integration.
- **Statistical Analysis**: Performs differential expression analysis with `limma`.
- **Visualization**: Generates high-quality plots, including:
  - Volcano plots with `EnhancedVolcano`
  - Heatmaps with `pheatmap`
  - Other diagnostic plots for data quality assessment

## Dependencies
The script uses the following R libraries:
- `GEOquery`
- `lumi`
- `limma`
- `EnhancedVolcano`
- `ggplot2`
- `pheatmap`
- And others (see the script for a full list)

## Dataset
The script is built to work with the GEO dataset **GSE113865**, but it can be adapted to other Illumina datasets with minimal modifications.

## Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/sananooor/illumina-data-analysis.git


Open the R script Finalscript-illumina.R in your R IDE (RStudio recommended).

   Ensure all required libraries are installed: R
   if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(c(
    "GEOquery", "lumi", "limma", "EnhancedVolcano", 
    "ggplot2", "pheatmap", "oligo", "oligoClasses",
    "illuminaHumanv3.db"))
      


Modify the script if necessary to accommodate other datasets.

Run the script to analyze the data and generate visualizations.

Output
Normalized and transformed data tables.

Diagnostic plots (e.g., PCA, boxplots).

Differential expression results.

Visualization files such as heatmaps and volcano plots.

## License
This project is licensed under the MIT License. See the LICENSE file for details.


