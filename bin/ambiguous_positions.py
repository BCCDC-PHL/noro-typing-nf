#!/usr/bin/env python3
import os, sys
import argparse
import pandas as pd 
import re
from collections import defaultdict

def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('pileup', help='Pileup file extracted from BAM alignment')
	parser.add_argument('--min_freq', default=0.2, help='Minimum frequency of ambiguous positions')
	parser.add_argument('--max_freq', default=0.7, help='Maximum frequency of ambiguous positions')
	parser.add_argument('-o', '--output', help='Pileup output of ambiguous positions only')
	return parser 

def find_ambiguous_single(pileup_df, min_freq, max_freq):
	
	frequencies = pileup_df[['count_A','count_C','count_G','count_T']].apply(lambda x: x / pileup_df['depth'])

	ambi_df = pileup_df.loc[((frequencies > min_freq) & (frequencies < max_freq)).any(axis=1)]

	return ambi_df


def main():
	parser = init_parser()
	args = parser.parse_args()
	pileup_df = pd.read_csv(args.pileup, sep='\t')

	ambi_df = find_ambiguous_single(pileup_df, args.min_freq, args.max_freq)
	
	if not args.output:
		outpath = re.split("_|\.", os.path.basename(args.pileup))[0]+ ".ambi.pileup.tsv"
	else:
		outpath = args.output
	ambi_df.to_csv(outpath, sep='\t', index=False)



if __name__ == '__main__':
	main()