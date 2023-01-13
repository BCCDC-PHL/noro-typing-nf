#!/usr/bin/env python3
#%%
import os 
import sys
import pandas as pd 
import argparse

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('blastn', help='BlastN results file in TSV format')
	parser.add_argument('-m', '--metric', default="bitscore", help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-r', '--ref_scores', required='bsr' in sys.argv, default=None, help='Raw self-blast scores of reference sequences used for Blast Score Ratio')
	parser.add_argument('-o','--output', default=None, help='Filtered results output')
	parser.add_argument('-c','--min_cov', default=25, type=int, help='Minimum coverage of database reference sequence by contig (percentage, default = 25)')
	parser.add_argument('-i','--min_id', default=90, type=int, help='Minimum nucleotide sequence identity between database reference sequence and contig (percentage, default = 90)')
	return parser

def parse_blast(filepath):
	cols = 'qseqid sseqid pident qlen slen bitscore rawscore'.split(' ')
	return pd.read_csv(filepath, sep='\t', names=cols)


def add_blast_score_ratio(blast_df, ref_score_path):
	ref_scores = pd.read_csv(ref_score_path, sep='\t')
	colnames = ref_scores.columns
	blast_df = blast_df.merge(ref_scores, left_on='sseqid', right_on=colnames[0]).drop(colnames[0],axis=1)
	blast_df['bsr'] = blast_df['rawscore'] / blast_df[colnames[1]]
	return blast_df

#%%
def filter_alignments(blast_results, score_column, min_cov, min_id):
	'''Find best contig for each genome segment. Returns datasheet with best contigs.'''
	print('Filtering alignments...')
	# Annotate alignments with segment and subtype
	# blast_results['segment'] = blast_results.apply(lambda row: row['sseqid'].split('|')[2], axis=1)
	blast_results[['genotype','strain']] = blast_results['sseqid'].str.split("|",expand=True)[[1,2]]
	# Discard alignments below minimum identity threshold
	blast_results = blast_results[blast_results['pident']>=min_id]
	# Keep only best alignment for each contig (by choice of scoring metric, i.e. bitscore or score)
	best_scores = blast_results[['qseqid', score_column]].groupby('qseqid').max().reset_index()
	blast_results = pd.merge(blast_results, best_scores, on=['qseqid', score_column])
	
	# ensure that final results only have a single genotype
	subtype_counts = blast_results[['qseqid', 'genotype']].drop_duplicates()
	subtype_counts = subtype_counts.groupby('qseqid').size().reset_index()
	subtype_counts = subtype_counts[subtype_counts[0]==1][['qseqid']]
	blast_results = pd.merge(blast_results, subtype_counts, on='qseqid')
	# Keep only alignments between contigs and ref seqs with median segment length
	median_slen = blast_results[['qseqid', 'slen']].groupby('qseqid').quantile(0.5, interpolation='higher').reset_index()
	blast_results = pd.merge(blast_results, median_slen, on=['qseqid', 'slen'])
	# Discard contigs that do not provide minimum coverage of a segment
	# Calculates the percent covered by the query relative to the subject
	blast_results = blast_results[blast_results['qlen'] * 100 / blast_results['slen'] >= min_cov]
	# De-duplicate sheet 
	blast_results = blast_results.drop_duplicates()
	return blast_results


def main():
	parser = get_parser()
	args = parser.parse_args()

	blast_df = parse_blast(args.blastn)

	if args.ref_scores:
		print("Computing blast score ratios...")
		blast_df = add_blast_score_ratio(blast_df, args.ref_scores)
		
	blast_df.to_csv(args.blastn.split(".tsv")[0] + '_bsr.tsv', index=False)

	df = filter_alignments(blast_df, args.metric, args.min_cov, args.min_id)

	if not args.blastn.endswith("tsv"):
		print("WARNING: Are blastn results in TSV format?")

	if not args.output:
		args.output = args.blastn.split(".tsv")[0] + '_filter.tsv'

	if df.shape[0] > 0:
		df.to_csv(args.output, index=False, sep='\t')
	else:
		print('WARNING: No valid contigs found.')
		df.to_csv(args.output, index=False, sep='\t')

if __name__ == '__main__':
	main()
# %%
