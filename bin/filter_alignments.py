#!/usr/bin/env python3
#%%
import os 
import sys
import pandas as pd 
import argparse

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('blastn', help='BlastN results file in TSV format')
	parser.add_argument('-o','--output', default='missing', help='Filtered results output')
	parser.add_argument('-c','--min_cov', default=25, type=int, help='Minimum coverage of database reference sequence by contig (percentage, default = 25)')
	parser.add_argument('-i','--min_id', default=90, type=int, help='Minimum nucleotide sequence identity between database reference sequence and contig (percentage, default = 90)')
	return parser

#%%
def filter_alignments(blast, min_cov, min_id):
	'''Find best contig for each genome segment. Returns datasheet with best contigs.'''
	print('Filtering alignments...')
	cols = 'qseqid sseqid pident bitscore qlen slen'.split(' ')
	blast_results = pd.read_csv(blast, sep='\t', names=cols)
	# Annotate alignments with segment and subtype
	# blast_results['segment'] = blast_results.apply(lambda row: row['sseqid'].split('|')[2], axis=1)
	blast_results['genotype'] = blast_results.apply(lambda row: row['sseqid'].split('|')[1], axis=1)
	blast_results['strain'] = blast_results.apply(lambda row: row['sseqid'].split('|')[2], axis=1)
	# Discard alignments below minimum identity threshold
	blast_results = blast_results[blast_results['pident']>=min_id]
	# Keep only best alingments for each contig (by bitscore)
	best_bitscores = blast_results[['qseqid', 'bitscore']].groupby('qseqid').max().reset_index()
	blast_results = pd.merge(blast_results, best_bitscores, on=['qseqid', 'bitscore'])
	# Discard contigs whose best alignments are to multiple segments
	# segment_counts = blast_results[['qseqid']].drop_duplicates()
	# #segment_counts = blast_results[['qseqid', 'segment']].drop_duplicates()
	# segment_counts = segment_counts.groupby('qseqid').size().reset_index()
	# segment_counts = segment_counts[segment_counts[0]==1][['qseqid']]
	# blast_results = pd.merge(blast_results, segment_counts, on='qseqid')
	# Discard contigs whose best alignments are to multiple subtypes
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

	df = filter_alignments(args.blastn, args.min_cov, args.min_id)

	if not args.blastn.endswith("tsv"):
		print("WARNING: Are blastn results in TSV format?")

	if args.output == 'missing':
		args.output = args.blastn.split(".tsv")[0] + '_filter.tsv'

	if df.shape[0] > 0:
		df.to_csv(args.output, index=False, sep='\t')
	else:
		print('WARNING: No valid contigs found.')
		df.to_csv(args.output, index=False, sep='\t')

if __name__ == '__main__':
	main()
# %%
