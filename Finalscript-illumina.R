# Load necessary libraries
library(GEOquery)
library(lumi)
library(Biobase)
library(oligo)
library(oligoClasses)
library(ggplot2)
library(gplots)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(geneplotter)
library(limma) 
library(dplyr)
library(stringr)
library(matrixStats) 
library(genefilter)
library(EnhancedVolcano)
library(illuminaHumanv3.db)

# Load the dataset from GEO
gset <- getGEO("GSE113865", GSEMatrix = TRUE)[[1]]

# Ensure the data dimensions match
exprs_data <- exprs(gset)
pheno_data <- pData(gset)
feature_data <- fData(gset)

# Add a new Phenotype column to pData based on source_name_ch1
pData(gset)$Phenotype <- ifelse(str_detect(pData(gset)$`source_name_ch1`, "tumor"), "tumor", "normal")

# Check dimensions
cat("Expression data dimensions: ", dim(exprs_data), "\n")
cat("Phenotype data dimensions: ", dim(pheno_data), "\n")
cat("Feature data dimensions: ", dim(feature_data), "\n")

# Ensure that the sample names match between expression data and phenotype data
if (!all(colnames(exprs_data) == rownames(pheno_data))) {
  cat("Sample names do not match between expression data and phenotype data.\n")
  stop("Sample names do not match between expression data and phenotype data")
}

# Ensure that the feature names match between expression data and feature data
if (!all(rownames(exprs_data) == rownames(feature_data))) {
  cat("Feature names do not match between expression data and feature data.\n")
  stop("Feature names do not match between expression data and feature data")
}

# Basic Quality Control
# Transposition and PCA construction of the log2 expression data
expression_data <- exprs_data
pca <- prcomp(t(expression_data), scale. = FALSE)

# Finding the percentage/ratio of the first 2 PCs
percentage <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
standard_deviation_ratio <- sqrt(percentage[2] / percentage[1])

pca_frame <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Phenotype = pData(gset)$Phenotype)
colnames(pca_frame)[2] <- "PC2"

# PCA Plot Visualization
pca_plot <- ggplot(data = pca_frame) + geom_point(mapping = aes(x = PC1, y = PC2, color = Phenotype))
pca_plot + labs(title = "PCA Plot of Log2 Transformed Expression Data") + 
  xlab(paste0("PC1, Exp: ", percentage[1], "%")) + 
  ylab(paste0("PC2, Exp: ", percentage[2], "%")) + 
  scale_shape_manual(values = c(4, 15)) + 
  scale_color_manual(values = c("red", "purple"))

ggsave(
  "pca_raw.png",
  plot = last_plot(),
  dpi = 900
)

# Box plot construction of Log2 Intensities of Raw Data
oligo::boxplot(gset, target = 'core', main = "Boxplot of Raw Expression Data")

# Summarization and Background Correction Without Normalization
raw_data_summarized <- lumiExpresso(gset, 
                                    bg.correct = TRUE, 
                                    bgcorrect.param = list(method = 'bgAdjust'), 
                                    variance.stabilize = FALSE,  # No scaling
                                    varianceStabilize.param = list(), 
                                    normalize = FALSE, 
                                    normalize.param = list(), 
                                    QC.evaluation = TRUE, 
                                    QC.param = list(), 
                                    verbose = TRUE)

row_medians <- Biobase::rowMedians(as.matrix(exprs(raw_data_summarized)))
rle <- sweep(Biobase::exprs(raw_data_summarized), 1, row_medians)
rle <- as.data.frame(rle)
rle_gathered <- gather(rle, sample_name, log2_expression_deviation)
rle_plot <- ggplot(rle_gathered, aes(sample_name, log2_expression_deviation))
rle_plot + geom_boxplot(outlier.shape = NA) + ylim(c(-2, 2)) + theme(axis.text.x = element_text(colour = "aquamarine4", angle = 60, size = 6.5, hjust = 1, face = "bold"))

# Summarization and Background Correction With Normalization
raw_data_normalized <- lumiExpresso(gset, 
                                    bg.correct = FALSE, 
                                    bgcorrect.param = list(method = 'bgAdjust'), 
                                    variance.stabilize = FALSE, 
                                    varianceStabilize.param = list(), 
                                    normalize = TRUE, 
                                    normalize.param = list(method = "quantile"), 
                                    QC.evaluation = TRUE, 
                                    QC.param = list(), 
                                    verbose = TRUE)

normalized_expression_data <- Biobase::exprs(raw_data_normalized)

# Box plot construction of Log2 Intensities of Raw Data
oligo::boxplot(normalized_expression_data, target = 'core', main = "Log2-Intensities of Normalized Data")

# Transposition and PCA construction of the normalised expression data
pca_normalized <- prcomp(t(normalized_expression_data), scale. = FALSE)

# Finding the percentage/ratio of the first 2 PCs
percentage_normalized <- round(100 * pca_normalized$sdev^2 / sum(pca_normalized$sdev^2), 1)
standard_deviation_ratio_normalized <- sqrt(percentage_normalized[2] / percentage_normalized[1])

pca_frame_normalized <- data.frame(PC1 = pca_normalized$x[,1], PC2 = pca_normalized$x[,2], Phenotype = pData(gset)$Phenotype)
colnames(pca_frame_normalized)[2] <- "PC2"

# PCA Visualisation
pca_plot_normalized <- ggplot(data = pca_frame_normalized) + geom_point(mapping = aes(x = PC1, y = PC2, color = Phenotype))
pca_plot_normalized + labs(title = "PCA Plot of Summarized, Calibrated and Normalized Data") + 
  xlab(paste0("PC1, Exp: ", percentage_normalized[1], "%")) + 
  ylab(paste0("PC2, Exp: ", percentage_normalized[2], "%")) + coord_fixed(ratio = standard_deviation_ratio_normalized) + scale_shape_manual(values = c(4, 15)) + scale_color_manual(values = c("red", "purple"))

ggsave(
  "pca_normalised.png",
  plot = last_plot(),
  dpi = 900
)

# Data Annotation for Heatmap Visualization
data_annotation <- data.frame(Phenotype = pData(gset)$Phenotype)
row.names(data_annotation) <- row.names(pData(raw_data_normalized))

distances <- as.matrix(dist(t(normalized_expression_data), method = "manhattan"))
rownames(distances) <- row.names(pData(raw_data_normalized))
hmcol <- rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(255))
colnames(distances) <- NULL
diag(distances) <- NA  # don't need any color in diagonal, because it's self-self
ann_colors <- list(Phenotype = c("normal" = "green", "tumor" = "red"))

map <- pheatmap(distances, annotation_row = data_annotation, annotation_colors = ann_colors, legend = TRUE, treeheight_row = 0, legend_breaks = c(min(distances, na.rm = TRUE), max(distances, na.rm = TRUE)), legend_labels = c("Small distance", "Large distance"), main = "Heatmap for Calibrated Samples (Normalized)")

# Intensity-based Filtering of Low-intensity Transcripts
threshold_set <- 2.8
transcript_medians <- rowMedians(Biobase::exprs(raw_data_normalized))
histo_transcript_medians <- hist(transcript_medians, 100, col = "lightblue", freq = FALSE,
                                 main = "Histogram of the Median Intensities Per Gene With Threshold",
                                 border = "black",
                                 xlab = "Median Intensities")
abline(v = threshold_set, col = "red", lwd = 2)

# Filtering out the Genes that are Above Threshold
no_of_samples <- table(pData(raw_data_normalized)$Phenotype)
samples_cutoff <- min(no_of_samples)
filtered_genes_data <- apply(Biobase::exprs(raw_data_normalized), 1,
                             function(x){
                               sum(x > threshold_set) >= samples_cutoff})
filtered_genes_table <- table(filtered_genes_data)
filtered_genes <- subset(raw_data_normalized, filtered_genes_data)

# Annotate the data with identifiers and names
library(illuminaHumanv4.db)

annotated_data <- AnnotationDbi::select(illuminaHumanv4.db,
                                        keys = rownames(filtered_genes),
                                        columns = c("SYMBOL", "GENENAME"),
                                        keytype = "PROBEID")

annotated_data <- subset(annotated_data, !is.na(SYMBOL))

# Removal of the PROBEIDS that match to multiple genes
grouped_genes <- dplyr::group_by(annotated_data, PROBEID)
grouped_genes_summarized <- dplyr::summarise(grouped_genes, total = n_distinct(SYMBOL))
filtered_annotated_genes <- filter(grouped_genes_summarized, total == 1)  # Keep probes with one gene
probe_statistics <- filtered_annotated_genes
dim(probe_statistics)

# Generation of ExpressionSet without the probe IDs that have multiple mappings
ids_to_exclude <- !(rownames(filtered_genes) %in% probe_statistics$PROBEID)
table(ids_to_exclude)
filtered_final_gene_set <- filtered_genes[!ids_to_exclude, ]
validObject(filtered_final_gene_set)

# Ensure that the filtered expression data has the correct dimensions
cat("Filtered expression data dimensions: ", dim(expression_data), "\n")

#Combining filtered data with symbols
fData(filtered_final_gene_set)$PROBEID <- rownames(fData(filtered_final_gene_set))
fData(filtered_final_gene_set) <- left_join(fData(filtered_final_gene_set), annotated_data)
rownames(fData(filtered_final_gene_set)) <-fData(filtered_final_gene_set)$PROBEID

write.csv(filtered_final_gene_set, "filtered_final_gene_set.csv")

# Differential Expression Analysis
# The new Phenotype column is used throughout
pData(filtered_final_gene_set)$phenotype <- ifelse(str_detect(Biobase::pData(filtered_final_gene_set)$Phenotype, "tumor"),"tumor", "normal")

#verify if the code added phenotype correctly
filtered_final_gene_set@phenoData@data[["Phenotype"]]

# Set up the design matrix
individual <- as.character(row.names(pData(filtered_final_gene_set)))
disease <- ifelse(str_detect(Biobase::pData(filtered_final_gene_set)$phenotype, "tumor"), "tumor", "normal")

design <- model.matrix(~ 0 + disease)
colnames(design)[1:2] <- c("normal", "tumor")
rownames(design) <- individual

design <- select(as.data.frame(design), `tumor`, `normal`)
design_sorted <- arrange(design, desc(`tumor`))

# Save the design matrix as a CSV file
write.csv(design_sorted, "design_matrix.csv")

# Fit the linear model
fit <- lmFit(filtered_final_gene_set, design)

# Set up contrasts (tumor vs normal)
contrast.matrix <- makeContrasts(tumor - `normal`, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract top differentially expressed genes
top_genes <- topTable(fit2, adjust = "fdr", sort.by = "P", number = Inf)

write.csv(top_genes, "differential_expression_results.csv")

# Separate upregulated and downregulated genes based on logFC and p-value cutoffs
upregulated_0_5 <- subset(top_genes, logFC > 0.5 & adj.P.Val < 0.05)
downregulated_0_5 <- subset(top_genes, logFC < -0.5 & adj.P.Val < 0.05)

upregulated_1 <- subset(top_genes, logFC > 1 & adj.P.Val < 0.05)
downregulated_1 <- subset(top_genes, logFC < -1 & adj.P.Val < 0.05)

# Save the separated genes as CSV files
write.csv(upregulated_0_5, "upregulated_0_5.csv")
write.csv(downregulated_0_5, "downregulated_0_5.csv")

write.csv(upregulated_1, "upregulated_1.csv")
write.csv(downregulated_1, "downregulated_1.csv")

# Visualization of Results
library(EnhancedVolcano)
# Create a volcano plot
EnhancedVolcano(top_genes, 
                lab = top_genes$SYMBOL,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "TNBC vs. Normal",
                gridlines.major = FALSE,
                gridlines.minor = FALSE)
