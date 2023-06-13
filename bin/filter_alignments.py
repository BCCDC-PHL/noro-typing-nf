#!/usr/bin/env python3
#%%
import os 
import sys
import pandas as pd 
import argparse
from tools import parse_fasta, write_fasta
import traceback

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('blastn', help='BlastN results file in TSV format')
	parser.add_argument('-m', '--metric', default="bitscore", help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-r', '--ref_scores', default=None, help='Raw self-blast scores of reference sequences used for Blast Score Ratio')
	parser.add_argument('-c','--contig_mode', action='store_true', help='Print out \
					best contig sequences instead of best references. References used by default.')
	parser.add_argument('-s', '--seqs', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('-t', '--header_pos_type', default=1, type=int, help='Zero-index position of the reference sequence type in the header(genotype or p-type)')
	parser.add_argument('-d', '--header_delim', default="|", help='Delimiter separating reference header fields.')
	parser.add_argument('-o','--tsv_out', default=None, help='Filtered results output')
	parser.add_argument('-O','--fasta_out', default=None, help='Filtered FASTA output')
	parser.add_argument('--min_cov', default=25, type=int, help='Minimum coverage of database reference sequence by contig (percentage, default = 25)')
	parser.add_argument('-i','--min_id', default=90, type=int, help='Minimum nucleotide sequence identity between database reference sequence and contig (percentage, default = 90)')
	return parser

def parse_blast(filepath, ref_scores, header_delim, header_pos_type):
	cols = 'qseqid sseqid pident qlen slen length bitscore rawscore'.split()
	blast_df = pd.read_csv(filepath, sep='\t', names=cols)

	try: 
		if blast_df.shape[0] == 0:
			raise ValueError("ERROR: BLAST input has 0 rows.")
		
		# split the header into genotype and strain columns 
		# blast_df[['genotype','strain']] = blast_df['sseqid'].str.split("|",expand=True)[[1,2]]
		# compute the coverage column
		# blast_df['coverage'] = blast_df['length'] * 100 / blast_df['slen']
		blast_df['prop_covered'] = blast_df['length'] * 100 / blast_df['slen']
		blast_df['type'] = blast_df['sseqid'].str.split(header_delim).str[header_pos_type]

		blast_df['cov'] = blast_df['qseqid'].str.split("_").str[-1].astype(float)
		#blast_df = blast_df.loc[blast_df['cov'] > 10]
		blast_df = blast_df.sort_values("cov", ascending=False)
		
		# BLAST SCORE RATIOS 
		if isinstance(ref_scores, pd.DataFrame):
			colnames = ref_scores.columns
			blast_df = blast_df.merge(ref_scores, left_on='sseqid', right_on=colnames[0]).drop(colnames[0],axis=1)
			blast_df['bsr'] = blast_df['rawscore'] * 100 / blast_df[colnames[1]]

		# add name column to the blast results 
		sample_name = os.path.basename(filepath).split("_")[0]
		blast_df.insert(0, 'sample_name', sample_name)
		return blast_df

	except Exception as e:
		print(traceback.format_exc())
		print("ERROR: Failed to filter BLAST outputs (error above). Creating empty output files.",file=sys.stderr)
		# final_cols = ['sample_name'] + cols + 'genotype strain prop_covered refscore bsr'.split()
		return pd.DataFrame()

#%%
def filter_alignments(blast_results, score_column, min_cov, min_id):
	'''
	Find best contig for each genome segment. Returns datasheet with best contigs.'''
	print('Filtering alignments...')
	# Annotate alignments with segment and subtype
	# blast_results['segment'] = blast_results.apply(lambda row: row['sseqid'].split('|')[2], axis=1)
	try: 
		blast_results = blast_results.drop_duplicates()

		score_cols = ['cov', 'bitscore', 'prop_covered']

		if score_column not in score_cols:
			score_cols += [score_column]

		scoring = blast_results[score_cols]
		scoring = (scoring - scoring.mean()) / scoring.std()
		scoring['composite'] = scoring.sum(axis=1)

		blast_results = pd.concat([blast_results, scoring[['composite']]], axis=1)

		# best_scores = blast_results[['qseqid', score_column]].groupby('qseqid')[score_column].nlargest(5).reset_index()[['qseqid',score_column]]
		# blast_results = blast_results.nlargest(5, score_column)

		idxmax = blast_results.groupby('type')['composite'].idxmax()
		filtered = blast_results.loc[idxmax].drop_duplicates('sseqid')
		filtered = filtered.sort_values('composite', ascending=False).reset_index(drop=True)

		if filtered.shape[0] > 3:
			filtered = filtered.loc[0:3]
		elif filtered.shape[0] < 3:
			rows_to_add = 3 - filtered.shape[0]
			na_df = pd.DataFrame([[None]*filtered.shape[1]] * rows_to_add, columns=filtered.columns)
			filtered = pd.concat([filtered, na_df])

		# blast_results = pd.merge(blast_results, best_scores, on=['qseqid', score_column])
		# blast_results = blast_results.sort_values(['qseqid', score_column],axis=0)
		
	except Exception as e:
		print(str(e))
		print("ERROR: Encountered an error while filtering BLAST results")

	return filtered

def write_best_sequence(blast_results, seq_path, fasta_out, mode):
	'''
	Looks up the single top hit in the DB FASTA file and writes them to their own FASTA file.
	'''

	if mode == 'reference':
		blast_field = 'sseqid'
	elif mode == 'contig':
		blast_field = 'qseqid'
	else:
		print("ERROR: Mode is neither reference or contig. Exiting. ", file=sys.stderr)
		sys.exit(1) 

	if blast_results.shape[0] == 0:
		print("ERROR: No BLAST results found when trying to write output sequence.", file=sys.stderr)
		sys.exit(1) 

	# filter to the single top hit 
	if blast_results.shape[0] > 1:
		blast_results = blast_results.iloc[0,:]

	# Open contigs FASTA and load seqs into dict (key=seq header)
	seqs = parse_fasta(seq_path)

	ids = [blast_results[blast_field]] if not isinstance(blast_results[blast_field], list) else blast_results[blast_field]
	# extract the top hit from the relevant sequence file 
	print( blast_results[blast_field])
	final_refs = {name : seqs[name] for name in ids}

	result = write_fasta(final_refs, fasta_out)

	return result

def main():
	parser = get_parser()
	args = parser.parse_args()

	# parse the optional reference scores file (needed for BLAST score ratio calculation)
	if args.ref_scores:
		ref_scores = pd.read_csv(args.ref_scores, sep='\t')
	else:
		ref_scores = None

	# parse the raw BLAST results 
	blast_df = parse_blast(args.blastn, ref_scores, args.header_delim, args.header_pos_type)
	print(blast_df)
	# exit if no BLAST results are found 
	if blast_df.shape[0] == 0:
		print("WARNING: No data found in blast input file. Exiting.", file=sys.stderr)
		sys.exit(1)
		# verbose_df = pd.DataFrame(columns=blast_df.columns)

	else:
		blast_df = filter_alignments(blast_df, args.metric, args.min_cov, args.min_id)

	if not args.tsv_out:
		args.tsv_out = args.blastn.split(".tsv")[0] + '_filter.tsv'

	# output blast results 
	blast_df.to_csv(args.tsv_out, index=False, sep='\t')
	# verbose_df.to_csv(args.tsv_out.split(".")[0] + "_full.tsv", index=False, sep='\t')

	# record which method is being used 
	mode = 'contig' if args.contig_mode else 'reference'

	# write the best FASTA outputs
	print(f"Running in {mode} mode.")
	result = write_best_sequence(blast_df, args.seqs, args.fasta_out, mode)
	print(f"Success: {result}")

if __name__ == '__main__':
	main()
# %%
