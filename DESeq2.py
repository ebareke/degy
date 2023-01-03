import argparse
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.pandas2ri as pandas2ri
from rpy2.robjects.packages import importr

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--count_matrix', required=True, help='Path to count matrix file')
  parser.add_argument('--sample_metadata', required=True, help='Path to sample metadata file')
  parser.add_argument('--log2_fold_change_cutoff', required=True, type=float, help='Log2 fold change cutoff for differential expression')
  parser.add_argument('--adjusted_p_value_cutoff', required=True, type=float, help='Adjusted p-value cutoff for differential expression')
  parser.add_argument('--output_prefix', required=True, help='Prefix for output files')
  return parser.parse_args()

def main():
  # parse command line arguments
  args = parse_args()
  count_matrix_file = args.count_matrix
  sample_metadata_file = args.sample_metadata
  log2_fold_change_cutoff = args.log2_fold_change_cutoff
  adjusted_p_value_cutoff = args.adjusted_p_value_cutoff
  output_prefix = args.output_prefix

  # read count matrix and sample metadata
  count_matrix = pd.read_csv(count_matrix_file, index_col=0)
  sample_metadata = pd.read_csv(sample_metadata_file, index_col=0)

  # check that count matrix and sample metadata have the same samples
  assert set(count_matrix.columns) == set(sample_metadata.index)

  # import DESeq2 R package
  DESeq2 = importr('DESeq2')

  # create DESeq2 object
  dds = DESeq2.DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_metadata, design='~ condition')

  # run differential expression analysis
  dds = DESeq2.DESeq(dds)

  # extract results as a data frame
  results = DESeq2.results(dds)
  results_df = pandas2ri.ri2py(results)

  # filter results by log2 fold change and adjusted p-value cutoffs
  significant_results_df = results_df[(np.abs(results_df['log2FoldChange']) > log2_fold_change_cutoff) & (results_df['padj'] < adjusted_p_value_cutoff)]
  significant_genes = significant_results_df.index.tolist()

  # check if there are any significant genes
  if len(significant_genes) > 0:
    # import gage R package
    gage = importr('gage')

    # perform enrichment analysis on significant genes
    gage_results = gage.gage(significant_genes, geneSetCollection='KEGG_2016', ref='h.all.v6.1.symbols', maxSize=1000)
   
        # extract gage results as a data frame
    gage_results_df = pandas2ri.ri2py(gage_results)

    # filter results by q-value cutoff
    q_value_cutoff = 0.05
    significant_gage_results_df = gage_results_df[gage_results_df['q.value'] < q_value_cutoff]

    # save significant gage results to file
    significant_gage_results_df.to_csv(f'{output_prefix}_significant_gage_results.tsv', sep='\t')

    # save enrichment plots
    robjects.r['pdf'](f'{output_prefix}_GO_biological_processes.pdf')
    gage.plot_gs(gage_results, geneSet='GO_Biological_Process', title='GO Biological Processes')
    robjects.r['dev.off']()

    robjects.r['pdf'](f'{output_prefix}_GO_molecular_functions.pdf')
    gage.plot_gs(gage_results, geneSet='GO_Molecular_Function', title='GO Molecular Functions')
    robjects.r['dev.off']()

    robjects.r['pdf'](f'{output_prefix}_KEGG_pathways.pdf')
    gage.plot_gs(gage_results, geneSet='KEGG_2016', title='KEGG Pathways')
    robjects.r['dev.off']()

if __name__ == '__main__':
  main()
