#!/usr/bin/env python3
import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re
import pandas as pd 

def init_parser():
	# create the top-level parser
	parser = argparse.ArgumentParser(prog='PROG')
	parser.add_argument('fasta', help='FASTA sequences to filter ')
	parser.add_argument('output', help='Output FASTA file')
	parser.add_argument('--header_delim', default='|', help='Delimiter character used in the reference FASTA header')
	parser.add_argument('--header_pos_accno', default=0, type=int, help='Zero-indexed position of accession number in FASTA header')
	parser.add_argument('--header_pos_type', default=1, type=int, help='Zero-indexed position of sequence type')
	parser.add_argument('-n','--min-prop-n', default=0.10, type=float, help='Minimum allowed proportion of Ns')
	parser.add_argument('-r','--rename', action='store_true', help='Rename sequence headers to a consistent schema used throughout the pipeline.')

	return parser 

def N_filter(seqs, min_prop_n):
	filtered = []
	removed = []
	for x in seqs:
		if x.count('N') / len(x) < min_prop_n:
			filtered.append(x)
		else:
			removed.append(x)
	return filtered, removed

def rename(seqs, header_delim, accno_pos, type_pos):

	for seq in seqs:

		fields = seq.id.split(header_delim)
		fields = [x.replace("_","-") for x in fields]
		extra = [x for n, x in enumerate(fields) if n not in {accno_pos, type_pos}]

		# reformat header into a consistent format used by the pipeline 
		# ACCNO|TYPE|...|
		newheader = "|".join([fields[accno_pos], fields[type_pos], *extra])
		seq.id = newheader
		seq.description = newheader
		seq.name = newheader

	return seqs

def remove_ambi(seqs):
	ambi = re.compile('[WRKYSMBDHV]')
	for s in seqs:
		s.seq = Seq(ambi.sub("N", str(s.seq)))
	return seqs

def main(args):
	seqs = list(SeqIO.parse(args.fasta, 'fasta'))

	# remove ambiguous bases 
	seqs = remove_ambi(seqs)

	# filter sequences based on completeness
	seqs, removed = N_filter(seqs, args.min_prop_n)
	
	# rename sequences to consistent schema if flag is provided
	if args.rename:
		seqs = rename(seqs, args.header_delim, args.header_pos_accno, args.header_pos_type)

	
	if len(removed) > 0:
		print("WARNING: The following sequences were removed:")
		for seq in removed:
			print(seq.id)
	
	SeqIO.write(seqs, args.outpath, 'fasta')

if __name__ == '__main__':
	parser = init_parser()
	args = parser.parse_args()

	main(args)