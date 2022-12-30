# Load required libraries
library(ggplot2)
library(plotly)
library(argparse)

# Set up argument parser
parser <- ArgumentParser(description='Parses arguments and generates a heatmap')
parser$add_argument('--count_file', type=character, required=TRUE, help='Tab-separated file with gene counts')
parser$add_argument('--test_column_idx', type=integer, required=TRUE, help='Indices of "test" condition columns')
parser$add_argument('--control_column_idx', type=integer, required=TRUE, help='Indices of "control" condition columns')
parser$add_argument('--out_prefix', type=character, required=TRUE, help='Prefix for output filename')

# Parse arguments
args <- parser$parse_args()

# Read in count file
counts <- read.delim(args$count_file, header=TRUE, row.names=1, sep='\t')

# Subset counts for test and control conditions
test_counts <- counts[, args$test_column_idx]
control_counts <- counts[, args$control_column_idx]

# Calculate fold changes between test and control conditions
fold_changes <- test_counts / control_counts

# Convert fold changes to log2 scale
log2_fold_changes <- log2(fold_changes)

# Generate heatmap using ggplot2
ggplot(data=log2_fold_changes, aes(x=rownames(log2_fold_changes), y=colnames(log2_fold_changes))) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low='yellow', high='red') +
  theme(axis.text.x=element_text(angle=90, hjust=1))

# Convert heatmap to plotly object
plotly_heatmap <- ggplotly()

# Save plotly object as png file
png(filename=paste0(args$out_prefix, '.png'))
plotly_heatmap
dev.off()
