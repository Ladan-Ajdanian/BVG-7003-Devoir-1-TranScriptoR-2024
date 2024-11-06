# Use-case
This script is designed for RNA-seq data analysis, specifically focusing on Cannabis samples. RNA-seq (RNA sequencing) is a powerful tool to measure gene expression by sequencing RNA molecules, and this analysis identifies genes that are differentially expressed between two conditions or groups (in this case, the XX and XY groups). Cannabis plants, like some animals, have a sexual chromosome system that determines their sex. Plants with XX chromosomes are female, producing female flowers that are particularly valuable in industrial and medicinal applications because they are rich in active compounds such as cannabinoids. On the other hand, plants with XY chromosomes are male, producing male flowers, which primarily play a role in pollen production and fertilization.


![image](https://github.com/user-attachments/assets/15ac7721-ead9-4841-af1d-c02d7ae6ca98)

The script performs several key steps in the analysis pipeline, including data normalization, statistical testing, visualization, and clustering. By analyzing these differences, researchers can gain insights into the biological mechanisms and functional differences between the experimental groups.

# Input Data
The script requires an RNA-seq dataset with the following format:

- Data Format: CSV file

- Gene IDs: Genes should be listed in rows.

- Samples: Each sample should be represented in columns.

- Additional Information: The first column should contain gene names, and the file should have numerical count data for each gene/sample combination.

* Example input file path: ~/BVG-Assignment/2_Data_RNASeq_Cannabis_Sex.csv

*Install Dependencies: Make sure that you have the necessary R packages installed: edgeR, ggplot2, pheatmap, and igraph.

```{r}
install.packages("edgeR")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("igraph")
```

  # Data Preprocessing
- Loading Data: The RNA-seq data, usually obtained from sequencing platforms, is loaded into R. The dataset is in the form of a matrix where rows represent genes and columns represent different samples. The first column contains gene identifiers.

- Defining Conditions: This step is essential for grouping the data based on experimental conditions. In this case, there are two groups: "XX" and "XY." These could represent different sex chromosomes in Cannabis plants (XX for female and XY for male, assuming a typical sex determination system). Defining these conditions allows for differential analysis to compare gene expression between these two groups.
  
```{r}
condition <- factor(c(rep("XX", 69), rep("XY", 69))) # Adjust according to your sample labels
group <- data.frame(condition)
```
The number 69 represents the number of samples in each experimental group being compared. Specifically, it indicates that there are 69 samples for group XX and 69 samples for group XY.

- DGEList Object Creation (edgeR): The DGEList object is a core object in edgeR (a package for RNA-seq analysis). It holds the count data (how many times each gene was observed in each sample) and the grouping information, which is crucial for subsequent analysis steps such as differential expression testing.

  ```{r}
dge <- DGEList(counts = data, group = condition)
```

- Filtering Low-Count Genes: In RNA-seq data, genes that have very low counts across all samples are often noise and not informative. Filtering out these genes improves the accuracy and power of the analysis by focusing on genes that are more likely to show meaningful biological variation.

```{r}
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
```

- Normalization (TMM): Normalization is performed to adjust for technical biases that can occur in RNA-seq data, such as differences in sequencing depth across samples. The TMM (Trimmed Mean of M-values) method is one approach that accounts for compositional differences between libraries and allows for more accurate comparison of gene expression levels.

```{r}
dge <- calcNormFactors(dge, method = "TMM")
normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)
```

- Differential Expression Analysis (edgeR): Differential expression analysis identifies genes that are significantly up- or down-regulated between the two experimental conditions. The statistical model used (GLM – Generalized Linear Model) assesses how gene expression differs between the XX and XY groups, adjusting for potential confounders like library size. A significance threshold (e.g., p-value < 0.05) is typically applied to identify genes of interest.

```{r}
dge <- estimateDisp(dge)
fit <- glmFit(dge, design = model.matrix(~condition))
lrt <- glmLRT(fit, coef = 2)  # Test for differences between "XX" and "XY"
```


  # Output & Instructions
  The script generates several outputs, which are useful for understanding gene expression patterns:

1- Significant Genes: A CSV file (significant_genes.csv) containing the gene list that is significantly differentially expressed between the two genders (XX vs. XY).

```{r}
significant_genes <- topTags(lrt, n = Inf, p.value = 0.05)$table
```

2- PCA (Principal Component Analysis) Plot: PCA is a dimensionality reduction technique that simplifies complex datasets by transforming them into principal components that explain the maximum variance in the data. It helps visualize patterns in the data and can reveal potential outliers, batch effects, or clustering patterns, such as differences between the XX and XY groups.

```{r}
norm_counts <- cpm(dge, log = TRUE)
```
```{r}
pca <- prcomp(t(norm_counts), scale. = TRUE)
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Condition = condition)
```
```{r}
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
geom_point(size = 3) + labs(title = "PCA Plot of Normalized Counts", x =
"PC1", y = "PC2") + theme_minimal() + scale_color_manual(values =
c("blue", "red"))
```

![001](https://github.com/user-attachments/assets/5f9119dc-eaf7-4d63-a868-0f858ec43e02)


3- Heatmap of Top Variable Genes: A heatmap is generated to visualize the expression levels of the top 500 most variable genes across the samples. Genes that show high variability between samples are often of great biological interest, as they might be involved in important regulatory processes or responses to experimental conditions.

```{r}
var_genes <- apply(norm_counts, 1, var)
top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:500]
top_var_data <- norm_counts[top_var_genes, ]
```
```{r}
annotation_data <- data.frame(Condition = condition)
rownames(annotation_data) <- colnames(top_var_data)
```
```{r}
if (nrow(top_var_data) > 1 && ncol(top_var_data) > 1) {
  pheatmap(
    top_var_data, scale = "row", show_rownames = FALSE,
    annotation_col = annotation_data,
    main = "Heatmap of Top Variable Genes",
    color = colorRampPalette(c("blue", "white", "red"))(50) # Custom color palette
  )
} else {
  message("Insufficient variable genes to generate a heatmap.")
}
```

![001](https://github.com/user-attachments/assets/22a05a5b-2116-4194-8c3b-0101914da05d)


4- Volcano Plot: A volcano plot combines both the fold-change (magnitude of expression difference between conditions) and statistical significance (p-value) of each gene. Genes that have both a large fold-change and high significance appear as "outliers" in the plot. This visualization helps to quickly identify genes that are strongly differentially expressed.

```{r}
significant_genes$logP <- -log10(significant_genes$PValue)
ggplot(significant_genes, aes(x = logFC, y = logP)) +
geom_point(aes(color = FDR < 0.05)) + scale_color_manual(values =
c("pink", "red")) + labs(title = "Volcano Plot", x = "Log Fold Change",
y = "-Log10 P-Value") + theme_minimal()
```

![001](https://github.com/user-attachments/assets/65a8dae8-3638-4862-83c1-64808a94c5d1)


5- Correlation Matrix: The correlation matrix shows how gene expression levels are correlated across samples. Highly correlated genes might be co-regulated or involved in the same biological pathway, while uncorrelated genes might represent independent processes. This matrix helps to assess the overall similarity or divergence between the different samples.

```{r}
correlation_matrix <- cor(normalized_counts)
pheatmap(correlation_matrix, main = "Correlation Matrix Heatmap", color
= colorRampPalette(c("blue", "white", "red"))(50))
```

![001](https://github.com/user-attachments/assets/82248fca-577a-480f-a238-b8a5f32f809e)


6-MA Plot: The MA plot visualizes the relationship between mean gene expression (M) and log fold-change (A) for each gene. This plot is useful for identifying genes with significant changes in expression, especially in the context of microarray or RNA-seq data, by visualizing both the intensity of expression and the direction of change.

```{r}
plotMD(lrt, column = 1, main = "MA Plot", xlab = "Average Log CPM", ylab = "Log Fold Change")
abline(h = c(-1, 1), col = "blue")
```

![001](https://github.com/user-attachments/assets/43d1daaf-ba7b-4ada-9160-1246b179ce18)


7- Gene Network Plot: A gene network visualizes the relationships or interactions between genes. In this case, significant genes identified in the differential expression analysis are used to construct the network. Genes that are closely related (either through physical interactions or co-expression) are connected in the network, providing insights into potential biological pathways or gene regulatory networks.

```{r}
gene_network <- graph_from_data_frame(d = significant_genes, directed = FALSE)
plot(gene_network, 
     vertex.size = 10, 
     vertex.label.cex = 0.7, 
     vertex.label.dist = 1,
     vertex.label.degree = -pi/2,
     main = "Gene Network")
```

![001](https://github.com/user-attachments/assets/2dc18d5c-ab72-4355-8a0e-0faede5d3352)


8- K-means Clustering: K-means clustering is an unsupervised machine learning algorithm used to group similar data points (or genes) based on their expression profiles. The number of clusters (K) is predefined, and the algorithm partitions the data into K clusters. This helps identify patterns or groups of genes that behave similarly across samples.

```{r}
set.seed(123)
kmeans_result <- kmeans(t(norm_counts), centers = 3)
```
```{r}
pca_data$Cluster <- as.factor(kmeans_result$cluster)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot with K-means Clustering", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))
```

![001](https://github.com/user-attachments/assets/0f26a88b-008b-43c6-baa1-f1c3ee4781b2)
