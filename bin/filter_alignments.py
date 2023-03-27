#!/usr/bin/env python3
#%%
import os 
import sys
import pandas as pd 
import argparse
from tools import parse_fasta, write_fasta
def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('blastn', help='BlastN results file in TSV format')
	parser.add_argument('-m', '--metric', default="bitscore", help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-r', '--ref_scores', required='bsr' in sys.argv, default=None, help='Raw self-blast scores of reference sequences used for Blast Score Ratio')
	parser.add_argument('-c','--contig_mode', action='store_true', help='Print out \
					best contig sequences instead of best references. References used by default.')
	parser.add_argument('-s', '--seqs', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('-o','--tsv_out', default=None, help='Filtered results output')
	parser.add_argument('-O','--fasta_out', default=None, help='Filtered FASTA output')
	parser.add_argument('--min_cov', default=25, type=int, help='Minimum coverage of database reference sequence by contig (percentage, default = 25)')
	parser.add_argument('-i','--min_id', default=90, type=int, help='Minimum nucleotide sequence identity between database reference sequence and contig (percentage, default = 90)')
	return parser

def parse_blast(filepath, ref_score_path):
	cols = 'qseqid sseqid pident qlen slen nident bitscore rawscore'.split()
	blast_df = pd.read_csv(filepath, sep='\t', names=cols)

	try: 
		if blast_df.shape[0] == 0:
			raise ValueError("ERROR: BLAST input has 0 rows.")
		
		# split the header into genotype and strain columns 
		blast_df[['genotype','strain']] = blast_df['sseqid'].str.split("|",expand=True)[[1,2]]
		# compute the coverage column
		# blast_df['coverage'] = blast_df['length'] * 100 / blast_df['slen']
		blast_df['prop_covered'] = blast_df['nident'] * 100 / blast_df['slen']
		
		# BLAST SCORE RATIOS 
		ref_scores = pd.read_csv(ref_score_path, sep='\t')
		colnames = ref_scores.columns
		blast_df = blast_df.merge(ref_scores, left_on='sseqid', right_on=colnames[0]).drop(colnames[0],axis=1)
		blast_df['bsr'] = blast_df['rawscore'] * 100 / blast_df[colnames[1]]

		# add name column to the blast results 
		sample_name = os.path.basename(filepath).split("_")[0]
		blast_df.insert(0, 'sample_name', sample_name)
		return blast_df

	except Exception as e:
		print(str(e))
		print("ERROR: Failed to filter BLAST outputs (error above). Creating empty output files.")
		final_cols = ['sample_name'] + cols + 'genotype strain prop_covered refscore bsr'.split()
		return pd.DataFrame(columns=final_cols)

#%%
def filter_alignments(blast_results, score_column, min_cov, min_id):
	'''
	Find best contig for each genome segment. Returns datasheet with best contigs.'''
	print('Filtering alignments...')
	# Annotate alignments with segment and subtype
	# blast_results['segment'] = blast_results.apply(lambda row: row['sseqid'].split('|')[2], axis=1)
	try: 
		# Discard alignments below minimum identity threshold
		# blast_results = blast_results[blast_results['pident']>=min_id]

		# blast_results = blast_results.loc[blast_results['prop_covered'] > 40]
		# Keep only best alignment for each contig (by choice of scoring metric, i.e. bitscore or score)

		blast_results = blast_results.drop_duplicates()
		best_scores = blast_results[['qseqid', score_column]].groupby('qseqid')[score_column].nlargest(5).reset_index()[['qseqid',score_column]]
		# blast_results = blast_results.nlargest(5, score_column)
		blast_results = pd.merge(blast_results, best_scores, on=['qseqid', score_column])
		blast_results = blast_results.sort_values(['qseqid', score_column],axis=0)

		verbose_results = blast_results.copy()

		# ensure that final results only have a single genotype
		# subtype_counts = blast_results[['qseqid', 'genotype']].drop_duplicates()
		# subtype_counts = subtype_counts.groupby('qseqid').size().reset_index()
		# subtype_counts = subtype_counts[subtype_counts[0]==1][['qseqid']]
		# blast_results = pd.merge(blast_results, subtype_counts, on='qseqid')
		# Keep only alignments between contigs and ref seqs with median segment length
		# median_slen = blast_results[['qseqid', 'slen']].groupby('qseqid').quantile(0.5, interpolation='higher').reset_index()
		# blast_results = pd.merge(blast_results, median_slen, on=['qseqid', 'slen'])
		# Discard contigs that do not provide minimum coverage of a segment
		# Calculates the percent covered by the query relative to the subject
		
		blast_results = blast_results.loc[[blast_results[score_column].idxmax()]]
		# De-duplicate sheet 
		

	except Exception as e:
		print(str(e))
		print("ERROR: Encountered an error while filtering BLAST results")


	return blast_results, verbose_results

def write_best_contigs(blast_results, contig_fasta, fasta_out):
	'''
	Looks up best contigs in contigs FASTA file and writes them to their own FASTA file.
	'''
	# De-duplicate rows from contigs with best alignments to multiple ref seqs 

	if blast_results.shape[0] == 0:
		open(fasta_out, 'w').close()
		return True

	if blast_results.shape[0] != 1:
		print("WARNING: Number of final contigs is not equal to 1.")

	# Open contigs FASTA and load seqs into dict (key=seq header)
	contig_seqs = parse_fasta(contig_fasta)

	final_contigs = {contig : contig_seqs[contig] for contig in blast_results['qseqid']}

	result = write_fasta(final_contigs, fasta_out)

	return result


def write_best_references(blast_results, ref_seqs_db, fasta_out):
	'''
	Looks up best ref seqs in ref seqs DB FASTA file and writes them to their own FASTA file.
	'''

	if blast_results.shape[0] == 0:
		open(fasta_out, 'w').close()
		return True

	if blast_results.shape[0] != 1:
		print("WARNING: Number of final contigs is not equal to 1.")

	# Open contigs FASTA and load seqs into dict (key=seq header)
	ref_seqs = parse_fasta(ref_seqs_db)

	final_refs = {ref_name : ref_seqs[ref_name] for ref_name in blast_results['sseqid']}

	result = write_fasta(final_refs, fasta_out)

	return result

def main():
	parser = get_parser()
	args = parser.parse_args()

	blast_df = parse_blast(args.blastn, args.ref_scores)

	if blast_df.shape[0] == 0:
		print("WARNING: No data found in blast input file. Exiting with empty outputs.")
		verbose_df = pd.DataFrame(columns=blast_df.columns)

	else:
		blast_df, verbose_df = filter_alignments(blast_df, args.metric, args.min_cov, args.min_id)


	if not args.tsv_out:
		args.tsv_out = args.blastn.split(".tsv")[0] + '_filter.tsv'

	blast_df.to_csv(args.tsv_out, index=False, sep='\t')
	verbose_df.to_csv(args.tsv_out.split(".")[0] + "_full.tsv", index=False, sep='\t')

	if args.contig_mode:
		print("Running in contig (assembly) mode.")
		result = write_best_contigs(blast_df, args.seqs, args.fasta_out)
	else:
		print("Running in reference (alignment) mode.")
		result = write_best_references(blast_df, args.seqs, args.fasta_out)

	print(f"Success: {result}")

if __name__ == '__main__':
	main()
# %%
