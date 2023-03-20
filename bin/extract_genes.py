#!/usr/bin/env python3
"""
Cutting out genes using pairwise alignment. 

This script uses pairwise alignments to extract genes from FASTA-formatted query sequences. 
There are two specific use cases described here: 

1. Given one query, multiple references, and associated gene positions, return the gene cut out of the query sequence.
2. Given multiple queries, one reference, and associated gene positions, return the genes of each query sequence in a single FASTA file. 

"""

import os, sys
from Bio import SeqIO, Align, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from copy import copy 
import argparse
import yaml
from tools import make_align_dict, translate_nuc, flex_translate, get_boundaries
import re

def init_parser():
	# parser = argparse.ArgumentParser()
	# parser.add_argument('-q', '--query', required=True, help='FASTA sequences to align to reference')
	# parser.add_argument('-r', '--ref', required=True, help='Reference sequence FASTA file')
	# parser.add_argument('-p', '--positions', required=True, help='YAML file containing gene positions for extraction')
	# parser.add_argument('-g', '--gene', default='INFER', help='Gene to extract. Only relevant for multi-query BLAST extraction.')
	# parser.add_argument('--debug', action='store_true')
	# parser.add_argument('--accno_pos', default=0, help='Zero indexed position of accession number in reference FASTA header')
	# parser.add_argument('--header_delim', default='|', help='Delimiter character used in the reference FASTA header')
	# parser.add_argument('--output', help='Output FASTA file containing genes')

	# create the top-level parser
	parser = argparse.ArgumentParser(prog='PROG')
	parser.add_argument('--debug', action='store_true')
	parser.add_argument('--accno_pos', default=0, help='Zero indexed position of accession number in FASTA header')
	parser.add_argument('--header_delim', default='|', help='Delimiter character used in the reference FASTA header')

	subparsers = parser.add_subparsers(help='sub-command help')
	# create the parser for the "a" command
	parser_a = subparsers.add_parser('database', help='One reference, multiple input queries. Extract genes from a BLAST DB of query inputs')
	parser_a.add_argument('-q', '--query', required=True, help='FASTA sequences to align to reference')
	parser_a.add_argument('-r', '--ref', required=True, help='Single reference sequence FASTA file. For norovirus, this is the G1 reference.')
	parser_a.add_argument('-p', '--positions', required=True, help='YAML file containing gene positions for extraction')
	parser_a.add_argument('-g', '--gene', required=True, help='Gene to extract.')
	parser_a.add_argument('-o', '--output', required=True, help='Output FASTA file containing genes')

	parser_a.set_defaults(func=main_database)

	# create the parser for the "b" command
	parser_b = subparsers.add_parser('sample', help='Single input query mode. Extract one gene from one sample.')
	parser_b.add_argument('-q', '--query', required=True, help='Single FASTA sequence to align to chosen reference')
	parser_b.add_argument('-r', '--ref', required=True, help='Multi FASTA file of reference sequences')
	parser_b.add_argument('-g', '--gene', required=True, help='Gene to extract. For norovirus, either vp1 or rdrp.')
	parser_b.add_argument('-o', '--output', required=True, help='Output FASTA file containing genes')
	parser_b.set_defaults(func=main_sample)

	return parser 

def init_aligner():
	aligner = Align.PairwiseAligner()
	aligner.mode = 'global'
	aligner.open_gap_score = -2.0
	aligner.extend_gap_score = -0.1
	aligner.target_end_gap_score = 0.0
	aligner.query_end_gap_score = 0.0
	
	return aligner

def extract_gene_by_position(aligner, ref, query, gene_pos, debug):

	# Perform pairwise alignment and get aligned sequences 
	align_iter = aligner.align(ref.seq, query.seq)
	ref_aligned, qry_aligned = next(align_iter)

	# Generate necessary alignment dictionaries to translate coordinates 
	align_dict, _ = make_align_dict(ref_aligned, qry_aligned)

	# STEP 1 - USE ALIGNMENT DICT TO TRANSLATE FROM REF POS TO ALIGNMENT POS
	align_start = align_dict[gene_pos[0]-1]
	align_end = align_dict[gene_pos[1]]  # JUSTIFIED BY SANITY CHECK 

	print(query.id)
	print(f"Reference Length: {len(ref)}")
	print(f"Query Length: {len(query)}")
	print(f"Start: {align_start}")
	print(f"End: {align_end}")

	# formulate the final gene extracted from the query 
	outseq = qry_aligned[align_start: align_end].replace("-","")
	print(f"Output Length: {len(outseq)} % 3 == {len(outseq) % 3}")

	if debug:
		print("REFERENCE:")
		print(ref_aligned)
		print("QUERY:")
		print(qry_aligned)
		print("REFERENCE:")
		print(ref_aligned[align_start: align_end+1])
		print("QUERY:")
		print(qry_aligned[align_start: align_end+1])
		print("REFERENCE:")
		print(translate_nuc(ref_aligned[align_start: align_end+1],0))
		print("QUERY:")
		print(translate_nuc(qry_aligned[align_start: align_end+1],0))
		print(f"Output Length: {len(outseq)} % 3 == {len(outseq) % 3}")

	outrecord = SeqIO.SeqRecord(
		Seq(outseq),
		id=query.id,
		name=query.name,
		description=query.description
	)

	return outrecord

def extract_gene_by_gaps(aligner, ref, query, debug):
	# Perform pairwise alignment and get aligned sequences 
	align_iter = aligner.align(ref.seq, query.seq)
	ref_aligned, qry_aligned = next(align_iter)

	print(query.id)
	align_start, align_end = get_boundaries(ref_aligned)

	print(f"Reference Length: {len(ref)}")
	print(f"Query Length: {len(query)}")
	print(f"Start: {align_start}")
	print(f"End: {align_end}")

	# find the  
	outseq = qry_aligned[align_start: align_end].replace("-","")
	print(f"Output Length: {len(outseq)} % 3 == {len(outseq) % 3}")

	if debug:
		print("REFERENCE:")
		print(ref_aligned)
		print("QUERY:")
		print(qry_aligned)
		print("REFERENCE:")
		print(ref_aligned[align_start: align_end])
		print("QUERY:")
		print(qry_aligned[align_start: align_end])
		print("REFERENCE:")
		print(translate_nuc(ref_aligned[align_start: align_end],0))
		print("QUERY:")
		print(translate_nuc(qry_aligned[align_start: align_end],0))

	outrecord = SeqIO.SeqRecord(
		Seq(outseq),
		id=query.id,
		name=query.name,
		description=query.description
	)

	return outrecord
def sanity_check(seq):

	test_seq = copy(seq)

	# remove excess bases that do not fit 
	if len(test_seq) % 3 != 0:
		extra_bases = len(test_seq) % 3
		test_seq.seq = test_seq.seq[:-extra_bases]


	aa_qry = test_seq.translate()

	if aa_qry.seq.count("*") > 1:
		print("WARNING: Translated nucleotide has multiple stop codons. Going to attempt flex translation.")
		print("REFERENCE: \n"+translate_nuc(str(seq.seq),0))
		print("QUERY: \n "+translate_nuc(str(seq.seq),0))
		
		_, stop_codon_count, frame_start = flex_translate(str(seq.seq))

		if stop_codon_count <= 1:
			print("FIXED.")
			seq.seq = seq.seq[frame_start:]

			return seq
			
		else:
			print("WARNING: Problematic sequence.")
			print("REFERENCE: \n"+translate_nuc(str(seq.seq),0))
			print("QUERY: \n"+translate_nuc(str(seq.seq),0))

			return seq

	return seq

def main_database(args):
	print("Extracting gene from database...")

	# load reference gene positions from the args.genes YAML file 
	with open(args.positions,'r') as infile:
		position_dict = yaml.safe_load(infile)

	aligner = init_aligner()

	# load the reference seq list 
	ref_seq = next(SeqIO.parse(args.ref, 'fasta'))

	# load the database (a collection of query sequences)
	qry_seqs = list(SeqIO.parse(args.query, 'fasta'))

	# extract the gene from all query sequences
	gene_records = list(map(lambda x : extract_gene_by_position(aligner, ref_seq, x, position_dict[args.gene], args.debug), qry_seqs))
	gene_records = [sanity_check(x) for x in gene_records]

	# check for failed cases 
	failed = sum((1 for x in gene_records if not isinstance(x, SeqRecord)))

	# report failed cases 
	if failed > 0:
		print(f"WARNING: {failed} problematic sequences have been generated. \
			Only writing good quality ones. Sequence output will be lower than input.")
		
		gene_records = [x for x in gene_records if isinstance(x, SeqRecord)]

	# produce an output filename if none is given 
	if not args.output:
		output = re.split("_|\.", os.path.basename(args.query))[0] + '.genes.fasta'
	else:
		output = args.output

	SeqIO.write(gene_records, output, 'fasta')

def main_sample(args):
	'''
	Main function for mode #1 (single query). Function used to extract a gene from a single query sequence 
	'''

	print("Extracting gene from sample...")
	aligner = init_aligner()

	# load a single query sequence
	qry_seq = next(SeqIO.parse(args.query, 'fasta'))

	# parse out the new header from the query sequence
	newheader = qry_seq.description.split()[-1]
	qry_seq.id = newheader
	qry_seq.name = newheader
	qry_seq.description = newheader

	# load in a collection of references (which are already cut down to gene of interest) (this is the BLAST DB FASTA)
	ref_seq_dict = SeqIO.to_dict(SeqIO.parse(args.ref,'fasta'))
	ref_seq_dict = {x.split(args.header_delim)[args.accno_pos]:y for x,y in ref_seq_dict.items()}

	# parse out the reference accession number from the query FASTA file 
	accnos = qry_seq.id.split("|")[2].split("_")

	# adaptive search for the correct reference sequence (avoids need for an extra parameter)
	ref_seq = None

	# this section can be expanded if more genes are desired
	if args.gene == 'vp1':
		ref_seq = ref_seq_dict[accnos[0]]
	elif args.gene == 'rdrp':
		ref_seq = ref_seq_dict[accnos[1]]
	else:
		print("ERROR: Not a valid gene entry.")

	# extract the gene and perform sanity checks for the right reading frame 
	gene_record = extract_gene_by_gaps(aligner, ref_seq, qry_seq, args.debug)
	gene_record = sanity_check(gene_record)

	# check whether the gene is properly extracted 
	if gene_record:
		if not args.output:
			output = re.split("_|\.", os.path.basename(args.query))[0]
		else:
			output = args.output

		SeqIO.write(gene_record, output, 'fasta')

	else:
		print("ERROR: Failed to extract gene from sequence", file=sys.stderr)
	
def main():
	parser = init_parser()
	args = parser.parse_args()

	# call the assigned default function chosen by the user input 
	args.func(args)

if __name__ == '__main__':
	main()