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
	parser.add_argument('--debug', action='store_true')
	parser.add_argument('--accno_pos', default=0, help='Zero indexed position of accession number in FASTA header')
	parser.add_argument('--header_delim', default='|', help='Delimiter character used in the reference FASTA header')

	subparsers = parser.add_subparsers(help='sub-command help')
	# create the parser for the "a" command
	parser_a = subparsers.add_parser('main', help='')
	parser_a.add_argument('fasta', help='FASTA sequences to filter ')
	parser_a.add_argument('output', help='Output FASTA file')
	parser_a.add_argument('--header_delim', default='|', help='Delimiter character used in the reference FASTA header')
	parser_a.add_argument('--header_pos_accno', default=0, type=int, help='Zero indexed position of accession number in FASTA header')
	parser_a.add_argument('--header_pos_type', default=1, type=int, help='Zero indexed position of sequence type')
	parser_a.add_argument('-n','--min-prop-n', default=0.10, type=float, help='Minimum allowed proportion of Ns')
	parser_a.add_argument('-r','--rename', action='store_true', help='Minimum allowed proportion of Ns')


	parser_a.set_defaults(func=main_standard)

	# create the parser for the "b" command
	parser_b = subparsers.add_parser('phylo', help='Filter sequences based on BLAST results')
	parser_b.add_argument('fasta', help='FASTA sequences to filter')
	parser_b.add_argument('blast', help='BLAST results used in filtering')
	parser_b.add_argument('output', help='Output FASTA file')
	parser_b.set_defaults(func=main_tree)

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
	extra = [x for n, x in enumerate(fields) if n not in {accno_pos, type_pos}]

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

def main_standard(args):
	seqs = list(SeqIO.parse(args.fasta, 'fasta'))

	seqs, removed = N_filter(seqs, args.min_prop_n)
	# seqs = remove_ambi(seqs)

	if args.rename:
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

def main_tree(args):
	blast_df = pd.read_csv(args.blast, sep='\t')
	
	seqs = list(SeqIO.parse(args.fasta, 'fasta'))
	seqs = {x.id.split("|")[0]: x for x in seqs}

if __name__ == '__main__':
	parser = init_parser()
	args = parser.parse_args()

	args.func(args)