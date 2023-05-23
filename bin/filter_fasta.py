#!/usr/bin/env python3
import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re

def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('fasta', help='FASTA sequences to filter ')
	parser.add_argument('output', help='Output FASTA file')
	parser.add_argument('--header_delim', default='|', help='Delimiter character used in the reference FASTA header')
	parser.add_argument('--header_pos_accno', default=0, type=int, help='Zero indexed position of accession number in FASTA header')
	parser.add_argument('--header_pos_type', default=1, type=int, help='Zero indexed position of sequence type')
	parser.add_argument('-n','--min-prop-n', default=0.35, type=float, help='Minimum allowed proportion of Ns')
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

def rename(seq, header_delim, accno_pos, type_pos):

	fields = seq.id.split(header_delim)
	fields = [x.replace("_","-") for x in fields]
	extra = [x for n, x in enumerate(fields) if n not in [accno_pos, type_pos]]

	newheader = "|".join([fields[accno_pos], fields[type_pos], *extra])
	seq.id = newheader
	seq.description = newheader
	seq.name = newheader

	return seq

def remove_ambi(seqs):
	ambi = re.compile('[WRKYSMBDHV]')
	for s in seqs:
		s.seq = Seq(ambi.sub("N", str(s.seq)))
	return seqs

def main():
	parser = init_parser()
	args = parser.parse_args()
	seqs = list(SeqIO.parse(args.fasta, 'fasta'))

	seqs, removed = N_filter(seqs, args.min_prop_n)
	# seqs = remove_ambi(seqs)
	seqs = list(map(lambda x: rename(x, args.header_delim, args.header_pos_accno, args.header_pos_type), seqs))

	if not args.output:
		outpath = os.path.basename(args.fasta).split(".")[0] + '_filtered.fasta'
	else:
		outpath = args.output
	
	if len(removed) > 0:
		print("WARNING: The following sequences were removed:")
		for seq in removed:
			print(seq.id)
	
	SeqIO.write(seqs, outpath, 'fasta')

if __name__ == '__main__':
	main()