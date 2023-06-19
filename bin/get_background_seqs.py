#!/usr/bin/env python3
#%%
import os, sys
from glob import glob 
import argparse 
import re
from Bio import SeqIO
import random 
from datetime import datetime as dt
def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('path', help='')
	parser.add_argument('-f', '--infasta', required=True, help='Multifasta file of current run sequences.')
	parser.add_argument('-s', '--schema', default='results*/phylo/{gene}/align/*.multi.fasta', help='Globbing path schema (with wildcards) used to search for consensus sequences. Use {gene} to indicate the variable gene name.')
	parser.add_argument('-g', '--gene', default='full', help='Type of sequences: 1) "full" genome 2) "capsid" gene or 3) "polymerase" gene')
	parser.add_argument('-o', '--outfasta', default='background_seqs.fasta', help='Name of the output FASTA file containing the background sequences')
	parser.add_argument('-m', '--maxseqs', default=50, type=int, help='Name of the output FASTA file containing the background sequences')
	return parser

def get_create_time(path):
	t = os.stat(path).st_ctime
	return dt.fromtimestamp(t)

def read_sequences(path, target_set):

	return [x for x in SeqIO.parse(path, 'fasta') if x.id.endswith('sample') and x.id.split("|")[0] not in target_set]

def sample_sequences(path, n, target_set):
	print(path, n)
	if n == 0:
		return []

	seqs = read_sequences(path, target_set)
	print(seqs)
	if n > len(seqs):
		return seqs
	else:
		return random.sample(seqs, n) 

def main():
	parser = init_parser()
	args = parser.parse_args()

	sample_ids = {x.id.split("|")[0] for x in SeqIO.parse(args.infasta, 'fasta') if 'sample' in x.id}
	print(sample_ids)
	# parse the file path schema
	schema = re.sub("\{.+\}", args.gene+"*", args.schema)

	# find all consensus FASTA files matching the glob schema
	consensus_files = glob(os.path.join(args.path, schema))

	# count the number of sequences per file 
	counts = list(map(lambda x: len(read_sequences(x, sample_ids)), consensus_files))

	max_sequences = args.maxseqs
	min_sequences = 0

	print(sum(counts))
	print(max_sequences)

	# downsampling 
	if sum(counts) > max_sequences:
		print("Downsampling sequences...")
		proportion = max_sequences / sum(counts) 

		# calculate proportion sample sizes 
		sample_counts = [int(x * proportion) for x in counts]
		# sample_counts = [x if x >= min_sequences else min_sequences for x in sample_counts]

		print("Original: ", sum(counts))
		print("Sampled: ", sum(sample_counts))

		# down sample to appropriate counts per file 
		seqs = [seq for path, n in zip(consensus_files, sample_counts) for seq in sample_sequences(path, n, sample_ids)]
	
	# no sampling required
	else:
		# simply extract all sequences
		seqs = [seq for path in consensus_files for seq in SeqIO.parse(path, 'fasta') if "sample" in seq.id and seq.id.split('|')[0] not in sample_ids]

	# rename sequences for background formatting
	for s in seqs:
		newheader = s.id.replace("sample", "background|"+dt.today().strftime('%Y-%m-%d'))
		s.id = newheader
		s.name = newheader
		s.description = newheader

	# write outputs 
	SeqIO.write(seqs, args.outfasta, 'fasta')
#%%
if __name__ == '__main__':
	main()
	
