# MIT License
#
# Copyright (c) 2023 Eric BAREKE (eric.bareke@mcgill.ca), Emma Carlson (emma.carlson@mail.mcgill.ca), and Majewski Lab (McGill University)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# How to run it ?
# Rscript DEA.r -f /path/to/your/input/files -c condition1,condition2 -r reference_condition -fc 1.5 -p 0.05 -m 10
# Rscript DEA.r -f /home/ebareke/Desktop/T718 -c Sensitive,Resistant -r Sensitive -fc 1.5 -p 0.05 -m 10

# List of required packages
libraries <- c(
  "clusterProfiler",
  "enrichplot",
  "ggplot2",
  "DESeq2",
  "EnhancedVolcano",
  "grDevices",
  "RColorBrewer",
  "pheatmap",
  "ggrepel",
  "dplyr",
  "optparse",
  "here",
  "org.Hs.eg.db",
  "pathfindR"
)

# Function to check, install, and load packages
check_install_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    # Install the package if not installed
    if (package %in% BiocManager::installed()) {
      BiocManager::install(package)
    } else {
      install.packages(package)
    }
  }
  
  # Load the package
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

# Apply the function to each package in the list
lapply(libraries, check_install_load)

# Function to parse command line arguments
parse_args <- function() {
  option_list <- list(
    make_option(c("-f", "--filePath"), type = "character", help = "Path to the input files"),
    make_option(c("-c", "--conditions"), type = "character", help = "Conditions to compare (comma-separated)"),
    make_option(c("-r", "--refStr"), type = "character", help = "Reference condition"),
    make_option(c("-fc", "--FC_cutoff"), type = "numeric", help = "Fold change cutoff", default = 1.5),
    make_option(c("-p", "--adjusted_p_value_cutoff"), type = "numeric", help = "Adjusted p-value cutoff", default = 0.05),
    make_option(c("-m", "--min_read_cutoff"), type = "numeric", help = "Minimum read count cutoff", default = 20)
  )
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
  return(opt)
}

# Function Definitions

# Function for saving plots
save_plot <- function(plot, filename) {
  # Define file_ext here using snake_case
  file_ext <- tools::file_ext(filename)

  ggsave(filename, plot = plot, path = outputPath, dpi = 600, device = 'png')
  ggsave(sub(paste0(".", file_ext), "_highres.", file_ext), plot = plot, path = outputPath, dpi = 1200, device = 'png')
}


# Perform inner joins
merge_and_rename <- function(data, subset_data, suffix) {
  merged_data <- merge(counts(dds, normalized = FALSE), subset_data, by = "row.names", all = FALSE)
  colnames(merged_data)[1] <- "GeneSymbol"
  colnames(merged_data)[-1] <- paste0(colnames(merged_data)[-1], suffix)
  return(merged_data)
}

# Write results to Excel-compatible files with raw counts and Gene Symbol
write_results <- function(data, filename) {
  write.table(data, file = here::here(outputPath, filename), sep = "\t", row.names = FALSE, quote = FALSE)
}

# Function to perform GO enrichment analysis
perform_GO_enrichment <- function(gene_list, ont_type, output_path, adjusted_p_value_cutoff) {
  GO_results <- enrichGO(
    gene = gene_list,
    universe = background,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont_type,
    pAdjustMethod = "fdr",
    minGSSize = 10,
    maxGSSize = 2000,
    pvalueCutoff = adjusted_p_value_cutoff,
    qvalueCutoff = adjusted_p_value_cutoff,
    readable = TRUE
  )

  # Plot and save charts...
  save_enrichment_plots(GO_results, "GO", ont_type, output_path)
}

perform_KEGG_enrichment <- function(gene_list, output_path, adjusted_p_value_cutoff) {
  KEGG_results <- enrichKEGG(
    gene = gene_list,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = adjusted_p_value_cutoff,
    pAdjustMethod = "BH",
    universe = background,
    minGSSize = 10,
    maxGSSize = 2000,
    qvalueCutoff = adjusted_p_value_cutoff,
    use_internal_data = FALSE
  )

  # Plot and save charts...
  save_enrichment_plots(KEGG_results, "KEGG", NULL, output_path)
}

# Save enrichment plots
save_enrichment_plots <- function(results, analysis_type, ont_type, output_path) {
  bar_plot <- barplot(results, showCategory = 20)
  ggsave(paste0("barPlot_", analysis_type, ifelse(!is.null(ont_type), paste0("_", ont_type), ""), "_", output_path),
         plot = bar_plot, dpi = 600, device = 'png')
  ggsave(paste0("barPlot_", analysis_type, ifelse(!is.null(ont_type), paste0("_", ont_type), ""), "_", output_path),
         plot = bar_plot, dpi = 600, device = 'pdf')

  dot_plot <- dotplot(results, showCategory = 20)
  ggsave(paste0("dotPlot_", analysis_type, ifelse(!is.null(ont_type), paste0("_", ont_type), ""), "_", output_path),
         plot = dot_plot, dpi = 600, device = 'png')
  ggsave(paste0("dotPlot_", analysis_type, ifelse(!is.null(ont_type), paste0("_", ont_type), ""), "_", output_path),
         plot = dot_plot, dpi = 600, device = 'pdf')
}


# Function to perform pathfindR analysis
perform_pathfindR_analysis <- function(input_data, output_path) {
  # Format data
  input_df <- data.frame(
    Gene.symbol = input_data$GeneSymbol,
    logFC = input_data$log2FoldChange,
    adj.P.Val = input_data$padj
  )

  # Define gene sets
  gene_sets <- c("KEGG", "GO-BP", "GO-MF")

  for (gene_set in gene_sets) {
    # Run pathfindR and save results
    output <- run_pathfindR(input_df, gene_sets = gene_set)
    clustered_output <- cluster_enriched_terms(output)

    # Create enrichment chart
    enrichment_chart(
      result_df = output,
      top_terms = 10
    )

    # Create clustered enrichment chart
    clustered_output <- enrichment_chart(clustered_output, plot_by_cluster = TRUE)
    
    # Save plots in different formats
    ggsave(paste0("clusterPlot_DN_", gene_set, "_", output_path), plot = clustered_output, dpi = 600, device = 'png')
    ggsave(paste0("clusterPlot_DN_", gene_set, "_", output_path), plot = clustered_output, dpi = 600, device = 'pdf')
  }
}


# Parse command line arguments
opt <- parse_args()

# Assign values to variables
filePath <- opt$filePath
conditions <- unlist(strsplit(opt$conditions, ","))
refStr <- opt$refStr
FC_cutoff <- opt$FC_cutoff
adjusted_p_value_cutoff <- opt$adjusted_p_value_cutoff
min_read_cutoff <- opt$min_read_cutoff

# Set input and output paths
inputPath <- file.path(filePath, "inputFiles")
outputPath <- file.path(filePath, "outputFiles")

# Create inputFiles folder if it doesn't exist
if (!dir.exists(inputPath)) {
  dir.create(inputPath, recursive = TRUE)
} else {
  cat("Input directory already exists. Skipping creation.\n")
}

# Create outputFiles folder if it doesn't exist
if (!dir.exists(outputPath)) {
  dir.create(outputPath, recursive = TRUE)
} else {
  cat("Output directory already exists. Skipping creation.\n")
}


# Move count matrix and sample information files to inputFiles folder
countMatrixFiles <- list.files(pattern = "*Counts.t*")
infoFiles <- list.files(pattern = "*Info.t*")

if (length(countMatrixFiles) > 0) {
  file.rename(countMatrixFiles, file.path(inputPath, countMatrixFiles))
}

if (length(infoFiles) > 0) {
  file.rename(infoFiles, file.path(inputPath, infoFiles))
}

countMatrixFile <- list.files(path = inputPath, pattern = ".*Counts.t.*")
sampleInfoFile <- list.files(path = inputPath, pattern = ".*Info.t.*")

countMatrix <- read.table(file.path(inputPath, countMatrixFile), header = TRUE, row.names = 1)
sampleInfo <- read.table(file.path(inputPath, sampleInfoFile), header = TRUE, row.names = 1)

# Remove the 'symbol' column (if any) from countMatrix
countMatrix <- countMatrix[, -which(names(countMatrix) == "symbol")]

# Filter rows based on minimum read count cutoff
countMatrix <- countMatrix[rowSums(countMatrix) >= min_read_cutoff, ]

# Subset countMatrix to include only samples present in sampleInfo
countMatrix <- countMatrix[colnames(countMatrix) %in% row.names(sampleInfo), ]

# DEseq - RNA seq analysis (Note: may need to change design element depending on comparison)
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = sampleInfo, design = ~ Condition)
dds$Condition <- relevel(dds$Condition, ref = refStr)

# Perform DESeq analysis
dds <- DESeq(dds)

# Access the normalized counts from DESeq2
normalized_counts <- counts(dds, normalized = TRUE)

# Calculate sample correlation matrix
sample_cor_matrix <- cor(normalized_counts)


# Create a heatmap using pheatmap
Samples <- pheatmap(sample_cor_matrix, 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         color = colorRampPalette(c("blue", "white", "red"))(20),
         main = "Sample Correlation Heatmap")

save_plot(Samples, "correlationPlot.png")
save_plot(Samples, "correlationPlot.pdf")


# Extract DESeq results
res <- results(dds)

# Filter significant DEGs based on adjusted p-value and fold change
resSig <- res[!is.na(res$padj) & (res$padj < adjusted_p_value_cutoff) & (abs(res$log2FoldChange) > log2(FC_cutoff)), ]

# Separate upregulated and downregulated genes
upregulatedResSig <- resSig[resSig$log2FoldChange > log2(FC_cutoff), ]
downregulatedResSig <- resSig[resSig$log2FoldChange < -log2(FC_cutoff), ]

merged_upregulatedResSig <- merge_and_rename(counts(dds, normalized = FALSE), upregulatedResSig, "_up")
merged_downregulatedResSig <- merge_and_rename(counts(dds, normalized = FALSE), downregulatedResSig, "_down")
merged_resSig <- merge_and_rename(counts(dds, normalized = FALSE), resSig, "_all")
merged_res <- merge_and_rename(counts(dds, normalized = FALSE), res, "_all")

# Write the merged dataframes to files
write_results(merged_upregulatedResSig, "upregulated_results.txt")
write_results(merged_downregulatedResSig, "downregulated_results.txt")
write_results(merged_resSig, "all_results.txt")
write_results(merged_res, "all_results_raw.txt")


# PCA
rld <- rlogTransformation(dds)
rv <- rowVars(assay(rld))
vsd <- varianceStabilizingTransformation(dds)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pc <- prcomp(t(assay(vsd)[select, ]))
condition <- sampleInfo$Condition
scores <- data.frame(pc$x, condition)
samples <- rownames(sampleInfo)
explained_variance <- summary(pc)$importance[2, ] * 100

AllSamples <- ggplot(scores, aes(x = PC1, y = PC2, col = factor(condition))) +
  xlab(paste0("PC1 Variance Explained: ", explained_variance[1], "%")) +
  ylab(paste0("PC2 Variance Explained: ", explained_variance[2], "%")) +
  geom_point(size = 5) +
  geom_label_repel(aes(label = samples), show.legend = FALSE, max.overlaps = 50) +
  ggtitle("Principal Component Analysis") +
  scale_colour_brewer(name = " ", palette = "Set1") +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    legend.position = "right",
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.y = element_line(linewidth = 0.5, color = "black"),
    axis.line.x = element_line(linewidth = 0.5, color = "black"),
    axis.text = element_text(color = "black"),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),
    panel.grid.major.x = element_line(color = "grey", linetype = "dashed")
  )

save_plot(AllSamples, "pcaPlot.png")
save_plot(AllSamples, "pcaPlot.pdf")


# EnrichedHeatMap

# Create custom_colors dynamically
custom_colors <- list(
  Condition = setNames(
    c("darkgreen", "magenta"),
    c(refStr, setdiff(conditions, refStr))
  )
)

countMatrix <- data.matrix(countMatrix)
top_varied_genes <- rownames(countMatrix)[order(rowVars(countMatrix, useNames = TRUE), decreasing = TRUE)][1:500]
MostVariedGenes <- pheatmap(
  countMatrix[top_varied_genes, ],
  main = "Top 500 most variable genes across samples",
  annotation_col = sampleInfo,
  annotation_colors = custom_colors,
  scale = "row",
  show_rownames = FALSE,
  legend_title = "Treatment Time",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)

save_plot(MostVariedGenes, "heatmapPlot.png")
save_plot(MostVariedGenes, "heatmapPlot.pdf")


# Volcano plot of up and down regulated genes
keyvals <- ifelse(
  (res$log2FoldChange < -log2(FC_cutoff) & res$padj < adjusted_p_value_cutoff), 'blue',
  ifelse((res$log2FoldChange > log2(FC_cutoff) & res$padj < adjusted_p_value_cutoff), 'red',
         'white')
)
keyvals[is.na(keyvals)] <- 'white'
names(keyvals)[keyvals == 'red'] <- 'Up-regulated'
names(keyvals)[keyvals == 'white'] <- ''
names(keyvals)[keyvals == 'blue'] <- 'Down-regulated'

# Assuming res is a data frame or a tibble with columns log2FoldChange and padj

# Calculate the range of log2FoldChange and padj
x_range <- range(res$log2FoldChange, na.rm = TRUE)
y_range <- range(res$padj, na.rm = TRUE)

# Set xlim and ylim dynamically
DEG <- EnhancedVolcano(res,
                                lab = "",
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'Differentially Expressed Genes',
                                subtitle = paste('Upregulated Genes:', sum(keyvals == 'red'),
                                                 'Downregulated Genes:', sum(keyvals == 'blue')),
                                caption = bquote("FC cutoff = .(FC_cutoff); p-adj cutoff = .(adjusted_p_value_cutoff)"),
                                xlab = bquote(~Log[2]~ 'Fold Change'),
                                xlim = x_range,
                                ylim = y_range,
                                pCutoff = 0.1,
                                FCcutoff = 0.0,
                                pointSize = 2.0,
                                colCustom = keyvals,
                                legendPosition = 'bottom',
                                legendLabSize = 14,
                                legendIconSize = 4.0,
                                gridlines.major = TRUE,
                                gridlines.minor = FALSE,
                                border = 'partial',
                                borderWidth = 1.2,
                                borderColour = 'black', max.overlaps = Inf
)

volcanoPlot <- volcanoPlot + theme(plot.subtitle = element_text(hjust = 0, size = 12))

save_plot(DEG, "volcanoPlot.png")
save_plot(DEG, "volcanoPlot.pdf")


# Define input parameters as needed...

universe <- merged_res$GeneSymbol
background <- mapIds(org.Hs.eg.db, keys = universe, column = "ENTREZID", keytype = "SYMBOL")

background <- background[!is.na(background)]

listDN <- merged_downregulatedResSig$GeneSymbol 
listDN <- mapIds(org.Hs.eg.db, keys = listDN, column = "ENTREZID", keytype = "SYMBOL")
listDN <- listDN[!is.na(names(listDN))]

# Main analysis for downregulated genes
perform_GO_enrichment(listDN, "BP", outputPath)
perform_GO_enrichment(listDN, "MF", outputPath)
perform_GO_enrichment(listDN, "CC", outputPath)
perform_KEGG_enrichment(listDN, outputPath)
perform_pathfindR_analysis(merged_downregulatedResSig, outputPath)

# Define other input parameters as needed

listUP <- merged_upregulatedResSig$GeneSymbol 
listUP <- mapIds(org.Hs.eg.db, keys = listUP, column = "ENTREZID", keytype = "SYMBOL")
listUP <- listUP[!is.na(names(listUP))]

# Main analysis for upregulated genes
perform_GO_enrichment(listUP, "BP", outputPath)
perform_GO_enrichment(listUP, "MF", outputPath)
perform_GO_enrichment(listUP, "CC", outputPath)
perform_KEGG_enrichment(listUP, outputPath)
perform_pathfindR_analysis(merged_upregulatedResSig, outputPath)
