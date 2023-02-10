#!/usr/bin/env python3
#%%
import os 
import sys
import pandas as pd 
import argparse
from Bio import SeqIO

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--gblast', required=True, help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-p', '--pblast', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('-G', '--gfasta', required=True, help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-P', '--pfasta', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	# parser.add_argument('--contig', action='store_true', help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	# parser.add_argument('--sample_name', help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	# parser.add_argument('--accno_pos', default=0, help='Zero indexed position of accession number')
	# parser.add_argument('--header_delim', default='|', help='Zero indexed position of accession number')
	parser.add_argument('-o', '--outfasta', default='best_reference.fasta', help='Output file in FASTA format.')
	parser.add_argument('-O', '--outblast', default='best_reference.fasta', help='Output file in FASTA format.')
	return parser

def select_best(gblast, pblast):
	if gblast['bsr'][0] > pblast['bsr'][0]:
		return "gtype"
	else:
		return "ptype"

# def rename(seq, sample_name, g_ref_name, p_ref_name, contig):
# 	stype = 'contig' if contig else 'ref'
# 	newname = '_'.join([sample_name, g_ref_name, p_ref_name, stype])
# 	seq.id = newname
# 	seq.name = newname
# 	seq.description = ''
# 	return seq

def main():

	parser = get_parser()
	args = parser.parse_args()

	gtype_df = pd.read_csv(args.gblast, sep='\t')
	ptype_df = pd.read_csv(args.pblast, sep='\t')

	gfasta = next(SeqIO.parse(args.gfasta, 'fasta'))
	pfasta = next(SeqIO.parse(args.pfasta, 'fasta'))

	if gtype_df.shape[0] != 1 or ptype_df.shape[0] != 1:
		print("WARNING: Filtered blast results do not contain a single entry as expected.")
		print(f"Genotype results: {gtype_df.shape[0]}")
		print(f"Ptype results: {ptype_df.shape[0]}")

	# if args.sample_name: 
	# 	sample_name = args.sample_name
	# else:
	# 	sample_name = os.path.basename(args.gblast).split("_")[0]
	
	# g_accno = gtype_df['sseqid'][0].split(args.header_delim)[args.accno_pos]
	# p_accno = ptype_df['sseqid'][0].split(args.header_delim)[args.accno_pos]

	if select_best(gtype_df, ptype_df) == 'gtype':
		# gfasta = rename(gfasta, sample_name, g_accno, p_accno, args.contig)
		SeqIO.write(gfasta, args.outfasta, 'fasta')
		gtype_df.to_csv(args.outblast, sep='\t', index=False)

	else:
		# pfasta = rename(pfasta, sample_name, g_accno, p_accno, args.contig)
		SeqIO.write(pfasta, args.outfasta, 'fasta')
		ptype_df.to_csv(args.outblast, sep='\t', index=False)



if __name__ == '__main__':
	main()
# %%
