#!/usr/bin/env python3
#%%
import os 
import sys
import pandas as pd 
import argparse
from Bio import SeqIO
import traceback

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--gblast', required=True, help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-p', '--pblast', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('-G', '--gfasta', required=True, help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	parser.add_argument('-P', '--pfasta', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	# parser.add_argument('--contig', action='store_true', help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	# parser.add_argument('--sample_name', help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('--accno_pos', default=0, help='Zero indexed position of accession number')
	parser.add_argument('--type_pos', default=1, help='Zero indexed position of sequence type')
	parser.add_argument('--header_delim', default='|', help='Zero indexed position of accession number')
	parser.add_argument('-o', '--outfasta', default='best_reference.fasta', help='Output file in FASTA format.')
	parser.add_argument('-O', '--outblast', default='best_reference.fasta', help='Output file in FASTA format.')
	return parser

def select_best(gblast, pblast):
	if gblast['bsr'][0] > pblast['bsr'][0]:
		return "gtype"
	else:
		return "ptype"

def rename(seq, header):
	# stype = 'contig' if contig else 'ref'
	# newname = '_'.join([sample_name, g_ref_name, p_ref_name, stype])
	# seq.id = header
	# seq.name = header
	seq.description = header
	return seq

def main():

	parser = get_parser()
	args = parser.parse_args()

	gtype_df = pd.read_csv(args.gblast, sep='\t')
	ptype_df = pd.read_csv(args.pblast, sep='\t')

	try:
		gfasta = next(SeqIO.parse(args.gfasta, 'fasta'))
	except Exception:
		print(traceback.format_exc())
		gfasta = None

	try:
		pfasta = next(SeqIO.parse(args.pfasta, 'fasta'))
	except Exception:
		print(traceback.format_exc())
		pfasta = None

	# [gtype, ptype]
	output_types = ['NA','NA']
	output_accno = ['NA','NA']

	# ERROR -- both outputs failed 
	if (gtype_df.shape[0] != 1 or not gfasta) and (ptype_df.shape[0] != 1 or not pfasta):
		print('ERROR: Problems with BOTH gtype and ptype outputs')
		print(f"Genotype results: {gtype_df.shape[0]}")
		print(f"Ptype results: {ptype_df.shape[0]}")
		sys.exit(1)
	
	# valid ptype outputs; gtyping failed
	elif gtype_df.shape[0] != 1 or not gfasta:
		choice = 'ptype'
		pfields = [x.replace("_","-") for x in pfasta.id.split(args.header_delim)]
		output_types[1] = pfields[args.type_pos]
		output_accno[1] =  pfields[args.accno_pos]
		sample_name = ptype_df['sample_name'][0]


	# valid gtype outputs; ptyping failed
	elif ptype_df.shape[0] != 1 or not pfasta:
		choice = 'gtype'
		gfields = [x.replace("_","-") for x in gfasta.id.split(args.header_delim)]
		output_types[0] = gfields[args.type_pos]
		output_accno[0] =  gfields[args.accno_pos]
		sample_name = gtype_df['sample_name'][0]

	
	# both outputs are valid
	else:
		choice = select_best(gtype_df, ptype_df)
		gfields = [x.replace("_","-") for x in gfasta.id.split(args.header_delim)]
		pfields = [x.replace("_","-") for x in pfasta.id.split(args.header_delim)]
		sample_name = gtype_df['sample_name'][0]
		output_types = [gfields[args.type_pos],	pfields[args.type_pos]]
		output_accno = [gfields[args.accno_pos],pfields[args.accno_pos]]

	newheader = "|".join([
		sample_name, 
		"_".join(output_types),
		"_".join(output_accno),
		"g" if choice == 'gtype' else 'p'
	])

	if choice == 'gtype':
		gfasta = rename(gfasta, newheader)
		SeqIO.write(gfasta, args.outfasta, 'fasta')
		gtype_df.to_csv(args.outblast, sep='\t', index=False)

	else:
		pfasta = rename(pfasta, newheader)
		SeqIO.write(pfasta, args.outfasta, 'fasta')
		ptype_df.to_csv(args.outblast, sep='\t', index=False)



if __name__ == '__main__':
	main()
# %%
