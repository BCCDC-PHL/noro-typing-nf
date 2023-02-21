#!/usr/bin/env python3
import os, sys
from Bio import SeqIO
import argparse

def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('fasta', help='FASTA sequences to filter ')
	parser.add_argument('output', help='Output FASTA file')
	parser.add_argument('-n','--min-prop-n', default=0.35, type=int, help='Minimum allowed proportion of Ns')
	return parser 

def filter(seqs, min_prop_n):
	seqs = [x for x in seqs if x.count('N') / len(x) < min_prop_n]
	return seqs

def main():
	parser = init_parser()
	args = parser.parse_args()
	seqs = list(SeqIO.parse(args.fasta, 'fasta'))

	filtered = filter(seqs, args.min_prop_n)

	if not args.output:
		outpath = os.path.basename(args.fasta).split(".")[0] + '_filtered.fasta'
	else:
		outpath = args.output
	SeqIO.write(filtered, outpath, 'fasta')

if __name__ == '__main__':
	main()