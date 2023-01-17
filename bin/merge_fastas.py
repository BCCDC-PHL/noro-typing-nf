#!/usr/bin/env python
import os, sys
import argparse

#%%
def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('fasta1')
	parser.add_argument('fasta2')
	parser.add_argument('output')
	parser.add_argument('--accno_pos', default=0, help="Zero-indexed position \
		where accession number is found in header. Default is first position (0)")
	parser.add_argument('--delim', default="|", help="Header delimiter")
	return parser

def write_fasta(seqs, outpath):
	with open(outpath, 'w') as outfile:
		for header, seq in seqs.items():
			outfile.write(">"+header+'\n')
			outfile.write(seq+'\n')

#%%
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

#%%
def union_fastas(fasta1, fasta2, delim, accno_pos):
	seqs1 = {x.split(delim)[accno_pos]:y for x,y in parse_fasta(fasta1).items()}
	seqs2 = {x.split(delim)[accno_pos]:y for x,y in parse_fasta(fasta2).items()}

	# union = set(seqs1.keys()).union(set(seqs2.keys()))
	union = seqs1 | seqs2

	# print('Length FASTA 1: ', len(seqs1), file=sys.stderr)
	# print('Length FASTA 2: ', len(seqs2), file=sys.stderr)
	# print('Length Union: ', len(union), file=sys.stderr)

	return union

if __name__ == '__main__':
	parser = init_parser()
	args = parser.parse_args()

	seqs = union_fastas(args.fasta1, args.fasta2, args.delim, args.accno_pos)

	write_fasta(seqs, args.output)

