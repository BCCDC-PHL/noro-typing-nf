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
	parser.add_argument('-g', '--gblast', required=True, help='Filtered BLAST results containing a single best hit.')
	parser.add_argument('-p', '--pblast', required=True, help='Filtered BLAST results containing a single best hit.')
	parser.add_argument('-G', '--gfasta', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('-P', '--pfasta', required=True, help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	# parser.add_argument('--contig', action='store_true', help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	# parser.add_argument('--sample_name', help='Sequences file in FASTA format. Either the contig file or reference BLAST database in FASTA format.')
	parser.add_argument('--header_pos_accno', default=0, help='Zero indexed position of accession number')
	parser.add_argument('--header_pos_type', default=1, help='Zero indexed position of sequence type')
	parser.add_argument('--header_delim', default='|', help='Delimiter used in incoming FASTA headers')
	parser.add_argument('-o', '--outfasta', default='best_reference.fasta', help='Output file in FASTA format.')
	parser.add_argument('-O', '--outblast', default='best_blast.tsv', help='Output BLAST results in TSV format.')
	return parser

def select_best(gblast, pblast, metric_column):
	if gblast[metric_column][0] > pblast[metric_column][0]:
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

	header_types = ['NA','NA']
	header_accnos = ['NA','NA']

	try:
		gfasta = next(SeqIO.parse(args.gfasta, 'fasta'))
		header_types[0] = gfasta.id.split(args.header_delim)[args.header_pos_type]
		header_accnos[0] = gfasta.id.split(args.header_delim)[args.header_pos_accno]

	except Exception:
		print(traceback.format_exc())
		gfasta = None

	try:
		pfasta = next(SeqIO.parse(args.pfasta, 'fasta'))
		header_types[1] = pfasta.id.split(args.header_delim)[args.header_pos_type]
		header_accnos[1] = pfasta.id.split(args.header_delim)[args.header_pos_accno]

	except Exception:
		print(traceback.format_exc())
		pfasta = None

	# BEST CASE -- both inputs are valid
	if (gtype_df.shape[0] > 0 and gfasta) and (ptype_df.shape[0] > 0 and pfasta):
		choice = select_best(gtype_df, ptype_df)
	
	# valid ptype outputs; gtyping failed
	elif ptype_df.shape[0] == 1 and pfasta:
		choice = 'ptype'

	# valid gtype outputs; ptyping failed
	elif gtype_df.shape[0] == 1 and gfasta:
		choice = 'gtype'
	
	# both outputs are valid
	else:
		print('ERROR: Problems with BOTH gtype and ptype outputs')
		print(f"Genotype results: {gtype_df.shape[0]}")
		print(f"Ptype results: {ptype_df.shape[0]}")
		sys.exit(1)

	output_seqs = [x for x in (gfasta, pfasta) if x]

	# add the full header as the description of all valid output reference seqs
	full_header = "_".join(header_types) + "|" + "_".join(header_accnos) + "|sample"
	for seq in output_seqs:
		seq.description = full_header

	# this finds the cases where there are identical genotype + ptype references  
	# only write one because this can cause the aligner to fail 
	if len(output_seqs) > 1 and len(set((x.id.split(args.header_delim)[args.header_pos_accno] for x in output_seqs))) == 1: 
		SeqIO.write(output_seqs[0], args.outfasta, 'fasta')

	# otherwise, write the outputs 
	elif len(output_seqs) > 0:
		SeqIO.write(output_seqs, args.outfasta, 'fasta')
	else:
		print("ERROR: Neither gtype or ptype reference was valid sequence. No FASTA output generated.")
		sys.exit(1)
	
	if choice == 'gtype':
		gtype_df.to_csv(args.outblast, sep='\t', index=False)
	else:
		ptype_df.to_csv(args.outblast, sep='\t', index=False)


if __name__ == '__main__':
	main()
# %%
