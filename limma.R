# Load required libraries
library(dplyr)
library(limma)
library(ggplot2)
library(clusterProfiler)
library(factoextra)

# Read in read count data and metadata files
count_data <- read.csv("read_count_data.csv", row.names = 1)
metadata <- read.csv("metadata.csv")

# Merge read count data with metadata
count_data_metadata <- count_data %>% left_join(metadata, by = "sample_id")

# Create a factor variable for control vs test samples
count_data_metadata$condition <- factor(count_data_metadata$condition, levels = c("control", "test"))

# Use limma to perform differential expression analysis
design <- model.matrix(~0 + condition, data = count_data_metadata)
colnames(design) <- levels(count_data_metadata$condition)
fit <- lmFit(count_data, design)
contrast.matrix <- makeContrasts(test - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract DEG results
deg_results <- topTable(fit2, coef = "test-control", number = Inf)

# Write DEG results to file
write.csv(deg_results, "deg_results.csv")

# Generate volcano plot
volcano_plot_data <- data.frame(log2FoldChange = deg_results$logFC,
                                -log10(deg_results$adj.P.Val))
ggplot(volcano_plot_data, aes(x = log2FoldChange, y = -log10(adj.P.Val))) +
  geom_point() +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")

# Save volcano plot to file
ggsave("volcano_plot.png", width = 8, height = 6)

# Generate PCA plot
pca <- prcomp(count_data, center = TRUE, scale. = TRUE)
fviz_pca_var(pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

# Save PCA plot to file
ggsave("pca_plot.png", width = 8, height = 6)

# Generate heatmap
hm_data <- count_data[,order(metadata$condition)]
hm_col <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
heatmap.2(hm_data,
          Colv = NA,
          dendrogram = "row",
          Rowv = NA,
          trace = "none",
          col = hm_col,
          margins = c(5, 10),
          cexRow = 0.8,
          cexCol = 0.8,
          key = TRUE,
          keysize = 1,
          density.info = "none",
          symkey = TRUE,
          notext = TRUE)

# Save heatmap to file
ggsave("heatmap_plot.png", width = 8, height = 6)
