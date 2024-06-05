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
	parser.add_argument('-m', '--metric', default="bitscore", type=str, help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-r', '--self_blast_scores', default=None, help='Raw self-blast scores of reference sequences used for Blast Score Ratio')
	parser.add_argument('-c','--contig_mode', action='store_true', help='Print out \
					best contig sequences instead of best references. References used by default.')
	parser.add_argument('-s', '--seqs', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('-t', '--header_pos_type', default=1, type=int, help='Zero-index position of the reference sequence type in the header(genotype or p-type)')
	parser.add_argument('-d', '--header_delim', default="|", help='Delimiter separating reference header fields.')
	parser.add_argument('-o','--tsv_out', default=None, help='Filtered results output')
	parser.add_argument('-O','--fasta_out', default=None, help='Filtered FASTA output')
	parser.add_argument('--min_cov', default=25, type=int, help='Minimum coverage of database reference sequence by contig (percentage, default = 25)')
	parser.add_argument('-i','--min_id', default=90, type=int, help='Minimum nucleotide sequence identity between database reference sequence and contig (percentage, default = 90)')
	parser.add_argument('-n','--top_n_hits', default=1, type=int, help='Define how many top hits to include per sample. Default is 1.')
	return parser

def load_blast_format(varname='BLAST_OUTFMT'):
	DEFAULT_OUTFMT = 'qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore score'
	blast_outfmt = os.environ.get(varname)

	if not blast_outfmt:
		print(f"BLAST_OUTFMT environment variable not found. Using default column format: {DEFAULT_OUTFMT}")
		blast_outfmt = DEFAULT_OUTFMT
	
	return blast_outfmt

def parse_blast(filepath):
	cols = load_blast_format().split()
	blast_df = pd.read_csv(filepath, sep='\t', names=cols)

	blast_df.insert(0, 'sample_name', os.path.basename(filepath).split("_")[0])
	return blast_df

def add_blast_columns(blast_df, self_blast_scores, header_delim, header_pos_type):

	if blast_df.shape[0] == 0:
		raise ValueError("ERROR: BLAST input has 0 rows.")
	
	blast_df['prop_covered'] = blast_df['length'] * 100 / blast_df['slen']
	blast_df['prop_covered'] = blast_df['prop_covered'].map(lambda x : x if x <= 100 else 100)   # to prevent coverages over 100%

	blast_df['type'] = blast_df['sseqid'].str.split(header_delim).str[header_pos_type]

	blast_df = blast_df.sort_values("bitscore", ascending=False)
	
	# BLAST SCORE RATIOS 
	if isinstance(self_blast_scores, pd.DataFrame):
		colnames = self_blast_scores.columns
		blast_df = blast_df.merge(self_blast_scores, left_on='sseqid', right_on=colnames[0]).drop(colnames[0],axis=1)
		blast_df['bsr'] = blast_df['rawscore'] * 100 / blast_df[colnames[1]]

	return blast_df

#%%
def filter_alignments(blast_results, score_columns, min_cov, min_id, top_n_hits):
	'''
	Find best contig for each genome segment. Returns datasheet with best contigs.'''
	print('Filtering alignments...')

	blast_results = blast_results.drop_duplicates()

	if 'bsr' in blast_results:
		score_columns += ['bsr']

	# apply a min-max scaling to the scoring columns
	scoring = blast_results[score_columns]
	scoring = (scoring - scoring.min()) / (scoring.max() - scoring.min())
	scoring['composite'] = scoring.sum(axis=1)

	blast_results = pd.concat([blast_results, scoring[['composite']]], axis=1)

	# group by contig and take only the best result
	idxmax = blast_results.groupby(['qseqid'])['composite'].idxmax()
	filtered = blast_results.loc[idxmax]

	filtered = filtered.sort_values('composite', ascending=False).reset_index(drop=True)

	if filtered.shape[0] > top_n_hits:
		filtered = filtered.iloc[0:top_n_hits,:]

	return filtered

def write_best_sequence(blast_results, seq_path, fasta_out, contig_mode):
	'''
	Looks up the single top hit in the DB FASTA file and writes them to their own FASTA file.
	'''

	if contig_mode:
		blast_field = 'qseqid'
	else:
		blast_field = 'sseqid'

	if blast_results.shape[0] == 0:
		print("ERROR: No BLAST results found when trying to write output sequence.", file=sys.stderr)
		sys.exit(1) 

	# filter to the single top hit 
	if blast_results.shape[0] > 1:
		blast_results = blast_results.iloc[0,:]

	# Open contigs FASTA and load seqs into dict (key=seq header)
	seqs = parse_fasta(seq_path)

	name = blast_results[blast_field] if isinstance(blast_results[blast_field], str) else blast_results[blast_field].tolist()[0]
	# extract the top hit from the relevant sequence file 

	final_refs = {name : seqs[name]}

	result = write_fasta(final_refs, fasta_out)

	return result

def main():
	parser = get_parser()
	args = parser.parse_args()

	# parse the optional reference scores file (needed for BLAST score ratio calculation)
	if args.self_blast_scores:
		self_blast_scores = pd.read_csv(args.self_blast_scores, sep='\t')
	else:
		self_blast_scores = None

	# parse the BLAST file
	blast_df = parse_blast(args.blastn)
	
	# exit if no BLAST results are found 
	if blast_df.shape[0] == 0:
		print("WARNING: No data found in blast input file. Exiting.", file=sys.stderr)
		sys.exit(1)

	blast_df = add_blast_columns(blast_df, self_blast_scores, args.header_delim, args.header_pos_type)
	
	blast_df = filter_alignments(blast_df, args.metric.split(","), args.min_cov, args.min_id, args.top_n_hits)

	if not args.tsv_out:
		args.tsv_out = args.blastn.split(".tsv")[0] + '_filter.tsv'

	# output blast results 
	blast_df.to_csv(args.tsv_out, index=False, sep='\t')
	# verbose_df.to_csv(args.tsv_out.split(".")[0] + "_full.tsv", index=False, sep='\t')

	# record which method is being used 
	mode = 'contig' if args.contig_mode else 'reference'

	# write the best FASTA outputs
	print(f"Running in {mode} mode.")
	result = write_best_sequence(blast_df, args.seqs, args.fasta_out, args.contig_mode)
	print(f"Success: {result}")

if __name__ == '__main__':
	main()
# %%
