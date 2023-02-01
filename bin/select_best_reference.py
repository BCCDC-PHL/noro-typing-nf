#!/usr/bin/env python3
#%%
import os 
import sys
import pandas as pd 
import argparse

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--gblast', required=True, help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-p', '--pblast', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	return parser

def select_best(gblast, pblast):
	if gblast['bsr'][0] > pblast['bsr'][0]:
		return "gtype"
	else:
		return "ptype"

def main(args):
	gtype_df = pd.read_csv(args.gblast, sep='\t')
	ptype_df = pd.read_csv(args.pblast, sep='\t')

	if gtype_df.shape[0] != 1 or ptype_df.shape[0] != 1:
		print("WARNING: Filtered blast results do not contain a single entry as expected.")
		print(f"Genotype results: {gtype_df.shape[0]}")
		print(f"Ptype results: {ptype_df.shape[0]}")

	print(select_best(gtype_df, ptype_df))


if __name__ == '__main__':
	parser = get_parser()
	args = parser.parse_args()

	main(args)
# %%
