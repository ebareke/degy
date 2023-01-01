#!/usr/bin/env Rscript

# Read and parse passed arguments
args <- commandArgs(trailingOnly = TRUE)
count_matrix_file <- args[1]
sample_metadata_file <- args[2]
lfc_cutoff <- as.numeric(args[3])
adjp_cutoff <- as.numeric(args[4])

# Read count matrix and sample metadata
count_matrix <- read.delim(count_matrix_file, row.names = 1)
sample_metadata <- read.delim(sample_metadata_file, row.names = 1)

# Perform differential expression analysis using DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_metadata, design = ~ condition)
dds <- DESeq(dds)
res_deseq2 <- results(dds)

# Perform differential expression analysis using edgeR
library(edgeR)
y <- DGEList(counts = count_matrix)
y <- calcNormFactors(y)
design <- model.matrix(~ condition, data = sample_metadata)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
res_edger <- topTags(lrt, n = nrow(count_matrix))

# Perform differential expression analysis using limma+voom
library(limma)
library(voom)
v <- voom(count_matrix, sample_metadata, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
res_limma <- topTable(fit, n = nrow(count_matrix))

# Save results for each approach
write.table(res_deseq2, "deseq2_results.txt", row.names = TRUE)
write.table(res_edger, "edger_results.txt", row.names = TRUE)
write.table(res_limma, "limma_results.txt", row.names = TRUE)

# Save significant results for each approach
sig_deseq2 <- res_deseq2[abs(res_deseq2$log2FoldChange) > lfc_cutoff & res_deseq2$padj < adjp_cutoff,]
write.table(sig_deseq2, "deseq2_sig_results.txt", row.names = TRUE)
sig_edger <- res_edger[abs(res_edger$logFC) > lfc_cutoff & res_edger$FDR < adjp_cutoff,]
write.table(sig_edger, "edger_sig_results.txt", row.names = TRUE)
sig_limma <- res_limma[abs(res_limma$logFC) > lfc_cutoff & res_limma$adj.P.Val < adjp_cutoff,]
write.table(sig_limma, "limma_sig_results.txt", row.names = TRUE)

# Produce and save QC plots for DESeq2
library(ggplot2)
library(pcaMethods)
library(gplots)

# PCA plot
pca <- prcomp(t(count_matrix), center = TRUE, scale = TRUE)
pc1 <- pca$x[,1]
pc2 <- pca$x[,2]
plot_data <- data.frame(pc1 = pc1, pc2 = pc2, condition = sample_metadata$condition)
ggplot(plot_data, aes(x = pc1, y = pc2, color = condition)) +
  geom_point(size = 2) +
  xlab(paste0("PC1 (", round(pca$sdev[1] / sum(pca$sdev) * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(pca$sdev[2] / sum(pca$sdev) * 100, 1), "%)"))
ggsave("deseq2_pca.png", width = 6, height = 6)

# Heatmap of top 500 most-varied genes
top_varied <- rowVars(count_matrix)[order(rowVars(count_matrix), decreasing = TRUE)][1:500]
count_matrix_top_varied <- count_matrix[top_varied,]
heatmap.2(as.matrix(count_matrix_top_varied), col = bluered(75), scale = "row", trace = "none",
          dendrogram = "row", main = "Top 500 Most-Varied Genes")
dev.copy(png, "deseq2_heatmap.png")
dev.off()

# Produce and save summary plots for DESeq2
# Scatter plot
ggplot(sig_deseq2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 2) +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value")
ggsave("deseq2_scatter.png", width = 6, height = 6)

# MA plot
ggplot(sig_deseq2, aes(x = log2FoldChange, y = log10(baseMean))) +
  geom_point(size = 2) +
  xlab("log2 Fold Change") +
  ylab("log10 Mean Expression")
ggsave("deseq2_ma.png", width = 6, height = 6)

# Volcano plot
ggplot(sig_deseq2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 2, alpha = 0.6) +
  xlim(c(-15, 15)) +
  ylim(c(0, 50)) +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value")
ggsave("deseq2_volcano.png", width = 6, height = 6)

# Produce and save summary plots for edgeR
# Scatter plot
ggplot(sig_edger, aes(x = logFC, y = -log10(FDR))) +
  geom_point(size = 2) +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value")
ggsave("edger_scatter.png", width = 6, height = 6)

# MA plot
ggplot(sig_edger, aes(x = logFC, y = log10(AveExpr))) +
  geom_point(size = 2) +
  xlab("log2 Fold Change") +
  ylab("log10 Mean Expression")
ggsave("edger_ma.png", width = 6, height = 6)

# Volcano plot
ggplot(sig_edger, aes(x = logFC, y = -log10(FDR))) +
  geom_point(size = 2, alpha = 0.6) +
  xlim(c(-15, 15)) +
  ylim(c(0, 50)) +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value")
ggsave("edger_volcano.png", width = 6, height = 6)

# Produce and save summary plots for limma+voom
# Scatter plot
ggplot(sig_limma, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 2) +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value")
ggsave("limma_scatter.png", width = 6, height = 6)

# MA plot
ggplot(sig_limma, aes(x = logFC, y = log10(AveExpr))) +
  geom_point(size = 2) +
  xlab("log2 Fold Change") +
  ylab("log10 Mean Expression")
ggsave("limma_ma.png", width = 6, height = 6)

# Volcano plot
ggplot(sig_limma, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 2, alpha = 0.6) +
  xlim(c(-15, 15)) +
  ylim(c(0, 50)) +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value")
ggsave("limma_volcano.png", width = 6, height = 6

#################################### Human enrichments ####################################################################
## Perform enrichment analysis on significantly expressed genes from DESeq2 using gprofiler2
#library(gprofiler2)
#
## Convert gene identifiers to ENTREZID
#sig_deseq2$ENTREZID <- mapIds(org.Hs.egSYMBOL, keys=sig_deseq2$gene, column=2, keytype="SYMBOL", multiVals="first")
#
## Select significantly expressed genes with ENTREZID
#sig_deseq2_entrez <- sig_deseq2[!is.na(sig_deseq2$ENTREZID),]
#
## Create a list of ENTREZIDs
#entrez_list <- as.list(sig_deseq2_entrez$ENTREZID)
#
## Perform enrichment analysis for GO biological processes, GO molecular functions, and KEGG pathways
#enrichment <- gprofiler2(query = entrez_list, organism = "hsapiens", user_threshold = 0.05, 
#                        ontologies = c("GO_BP", "GO_MF", "KEGG_PATHWAY"))
#
## Filter for significant enrichments
#sig_enrichment <- enrichment[enrichment$adj.p.value < adjp_cutoff,]
#
## Save enrichment results
#write.table(sig_enrichment, file="deseq2_enrichment.txt", sep="\t", row.names=FALSE)
#
## Produce and save summary plots
#library(ggplot2)
#
## Barplot of enrichment results
#ggplot(sig_enrichment, aes(x=reorder(term, p.value), y=p.value)) +
#  geom_bar(stat="identity", fill="darkblue") +
#  xlab("Term") +
#  ylab("p-value") +
#  coord_flip()
#ggsave("deseq2_enrichment.png", width = 6, height = 6)
#
## Repeat for the other approaches using sig_edger_entrez and sig_limma_entrez


# Perform enrichment analysis on significantly expressed genes from DESeq2 using gprofiler2
library(gprofiler2)

# Convert gene identifiers to ENTREZID
sig_deseq2$ENTREZID <- mapIds(org.Mm.egSYMBOL, keys=sig_deseq2$gene, column=2, keytype="SYMBOL", multiVals="first")

# Select significantly expressed genes with ENTREZID
sig_deseq2_entrez <- sig_deseq2[!is.na(sig_deseq2$ENTREZID),]

# Create a list of ENTREZIDs
entrez_list <- as.list(sig_deseq2_entrez$ENTREZID)

# Perform enrichment analysis for GO biological processes, GO molecular functions, and KEGG pathways
enrichment <- gprofiler2(query = entrez_list, organism = "mmusculus", user_threshold = 0.05, 
                        ontologies = c("GO_BP", "GO_MF", "KEGG_PATHWAY"))

# Filter for significant enrichments
sig_enrichment <- enrichment[enrichment$adj.p.value < adjp_cutoff,]

# Save enrichment results
write.table(sig_enrichment, file="deseq2_enrichment.txt", sep="\t", row.names=FALSE)

# Produce and save summary plots
library(ggplot2)

# Barplot of enrichment results
ggplot(sig_enrichment, aes(x=reorder(term, p.value), y=p.value)) +
  geom_bar(stat="identity", fill="darkblue") +
  xlab("Term") +
  ylab("p-value") +
  coord_flip()
ggsave("deseq2_enrichment.png", width = 6, height = 6)

# Perform enrichment analysis on significantly expressed genes from edgeR using gprofiler2

# Convert gene identifiers to ENTREZID
sig_edger$ENTREZID <- mapIds(org.Mm.egSYMBOL, keys=sig_edger$gene, column=2, keytype="SYMBOL", multiVals="first")

# Select significantly expressed genes with ENTREZID
sig_edger_entrez <- sig_edger[!is.na(sig_edger$ENTREZID),]

# Create a list of ENTREZIDs
entrez_list <- as.list(sig_edger_entrez$ENTREZID)

# Perform enrichment analysis for GO biological processes, GO molecular functions, and KEGG pathways
enrichment <- gprofiler2(query = entrez_list, organism = "mmusculus", user_threshold = 0.05, 
                        ontologies = c("GO_BP", "GO_MF", "KEGG_PATHWAY"))

# Filter for significant enrichments
sig_enrichment <- enrichment[enrichment$adj.p.value < adjp_cutoff,]

# Save enrichment results
write.table(sig_enrichment, file="edger_enrichment.txt", sep="\t", row.names=FALSE)

# Produce and save summary plots

# Barplot of enrichment results
ggplot(sig_enrichment, aes(x=reorder(term, p.value), y=p.value)) +
  geom_bar(stat="identity", fill="darkblue") +
  xlab("Term") +
  ylab("p-value") +
  coord_flip()
ggsave("edger_enrichment.png", width = 6, height = 6)

# Perform enrichment analysis on significantly expressed genes from limma+voom using gprofiler2

# Convert gene identifiers to ENTREZID
sig_limma$ENTREZID <- mapIds(org.Mm.egSYMBOL, keys=sig_limma$gene, column=2, keytype="SYMBOL", multiVals="first")

# Select significantly expressed genes with ENTREZID
sig_limma_entrez <- sig_limma[!is.na(sig_limma$ENTREZID),]

# Create a list of ENTREZIDs
entrez_list <- as.list(sig_limma_entrez$ENTREZID)

# Perform enrichment analysis for GO biological processes, GO molecular functions, and KEGG pathways
enrichment <- gprofiler2(query = entrez_list, organism = "mmusculus", user_threshold = 0.05, 
                        ontologies = c("GO_BP", "GO_MF", "KEGG_PATHWAY"))

# Filter for significant enrichments
sig_enrichment <- enrichment[enrichment$adj.p.value < adjp_cutoff,]

# Save enrichment results
write.table(sig_enrichment, file="limma_enrichment.txt", sep="\t", row.names=FALSE)

# Produce and save summary plots

# Barplot of enrichment results
ggplot(sig_enrichment, aes(x=reorder(term, p.value), y=p.value)) +
  geom_bar(stat="identity", fill="darkblue") +
  xlab("Term") +
  ylab("p-value") +
  coord_flip()
ggsave("limma_enrichment.png", width = 6, height = 6)

