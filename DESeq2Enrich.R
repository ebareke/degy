## Load required packages
library(DESeq2)
library(clusterProfiler)
library(enrichplot)

## Parse arguments
args <- commandArgs(trailingOnly=TRUE)
log2_fold_change_cutoff <- as.numeric(args[1])
adjusted_p_value_cutoff <- as.numeric(args[2])
count_matrix_file <- args[3]
metadata_file <- args[4]
organism <- args[5]

## Read count matrix and metadata as data frames
count_matrix <- read.delim(count_matrix_file, sep="\t", header=TRUE)
metadata <- read.delim(metadata_file, sep="\t", header=TRUE)

## Convert count matrix to a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=metadata, design= ~ condition)

## Run differential expression analysis
dds <- DESeq(dds)

## Extract results and apply log2 fold change and adjusted p-value cutoffs
results <- results(dds, lfcThreshold = log2_fold_change_cutoff, pvalThreshold = adjusted_p_value_cutoff)

## Extract significant genes
significant_genes <- rownames(results)[which(results$padj <= adjusted_p_value_cutoff)]

## Perform GO and KEGG enrichment analyses
go_enrichment <- enrichGO(significant_genes, organism=organism)
kegg_enrichment <- enrichKEGG(significant_genes, organism=organism)

## Create dotplot and save to file
dotplot(go_enrichment, showCategory=5)
ggsave("go_dotplot.pdf", plot=last_plot())

## Create enrichment map and save to file
enrichment_map(go_enrichment)
ggsave("go_enrichment_map.pdf", plot=last_plot())

## Create ridgeplot and save to file
ridgeplot(kegg_enrichment)
ggsave("kegg_ridgeplot.pdf", plot=last_plot())
