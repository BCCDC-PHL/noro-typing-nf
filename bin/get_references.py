#!/usr/bin/env python3
#%%
import os 
import sys
import pandas as pd 
import argparse

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-b','--blast', required=True, help='BlastN results file in TSV format')
	parser.add_argument('-s', '--seqs', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('-c','--contig_mode', action='store_true', help='Print out \
						best contig sequences instead of best references. References used by default.')
	parser.add_argument('-o','--output', default='missing', help='Filtered results output')
	return parser

def parse_fasta(filepath):	
	seqs = {}
	with open(filepath, 'r') as handle:
		for line in handle.readlines():
			if line[0] == '>':
				header = line.strip().lstrip('>')
				seqs[header] = ''
			else:
				seqs[header] += line.strip()
	return seqs

def parse_blastn(filepath):
	blast_results = pd.read_csv(filepath, sep='\t')
	blast_results = blast_results.drop_duplicates()

	sample_name = os.path.basename(filepath).split('_')[0]
	return blast_results, sample_name


def write_best_contigs(blast_path, contig_fasta, output_file):
	'''
	Looks up best contigs in contigs FASTA file and writes them to their own FASTA file.
	'''
	# De-duplicate rows from contigs with best alignments to multiple ref seqs 
	blast_results, sample_name = parse_blastn(blast_path)

	# Open contigs FASTA and load seqs into dict (key=seq header)
	contig_seqs = parse_fasta(contig_fasta)

	try: 
		with open(output_file, 'w') as outfile:

			for name, row in blast_results.iterrows():
				ref_name = row['qseqid']
				header = f'>{sample_name}_seq_{str(name + 1)}' \
						f"|{row['genotype']}|{row['strain']}|{row['slen']}|contig"
				outfile.write(header + '\n')
				outfile.write(contig_seqs[ref_name] + '\n')
	except Exception as e:
		print(str(e))
		return False

	return True

def write_best_references(blast_path, ref_seqs_db, output_file):
	'''
	Looks up best ref seqs in ref seqs DB FASTA file and writes them to their own FASTA file.
	'''
	# De-duplicate rows from contigs with best alignments to multiple ref seqs 
	blast_results, sample_name = parse_blastn(blast_path)
	
	best_bitscores = blast_results[['genotype', 'bitscore']].groupby(['genotype']).max().reset_index()
	blast_results = pd.merge(blast_results, best_bitscores, on=['genotype','bitscore'])

	# Chose ref seqs with median length for each segment/subtype combination
	median_lengths = blast_results[['genotype', 'slen']].groupby(['genotype']).quantile(0.5, interpolation='higher').reset_index()
	blast_results = pd.merge(blast_results, median_lengths, on=['genotype', 'slen'])
	# Choose first alphabetical ref seq for each segment/subtype combination
	first_ref_seqs = blast_results[['genotype', 'sseqid']].groupby(['genotype']).min().reset_index()
	blast_results = pd.merge(blast_results, first_ref_seqs, on=['genotype', 'sseqid'])

	# Open contigs FASTA and load seqs into dict (key=seq header)
	ref_seqs = parse_fasta(ref_seqs_db)

	# Write best contigs to FASTA
	try: 
		with open(output_file, 'w') as outfile:

			for name, row in blast_results.iterrows():
				ref_name = row['sseqid']
				header = f'>{sample_name}_seq_{str(name + 1)}' \
						f"|{row['genotype']}|{row['strain']}|{row['slen']}|ref"
				outfile.write(header + '\n')
				outfile.write(ref_seqs[ref_name] + '\n')
	except Exception as e:
		print(str(e))
		return False
	
	return True

if __name__ == '__main__':
	parser = get_parser()
	args = parser.parse_args()

	if args.output == 'missing':
		args.output = args.blast.split("_")[0] + '.refs.fasta'
	
	if args.contig_mode:
		print("Running in contig (assembly) mode.")
		result = write_best_contigs(args.blast, args.seqs, args.output)
	else:
		print("Running in reference (alignment) mode.")
		result = write_best_references(args.blast, args.seqs, args.output)
		

	if result:
		print(f"COMPLETE. Best sequence(s) written to {args.output}")
	else:
		print("ERROR: Failed to write all contigs.")
# %%
