#!/usr/bin/env python3
import os, sys
from Bio import SeqIO, Align, AlignIO
from Bio.Seq import Seq
import argparse
import yaml
from tools import make_align_dict, translate_nuc, flex_translate
from functools import partial
import multiprocessing as mp

def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('query', help='FASTA sequences to align to reference ')
	# parser.add_argument('g_ref', help='Genotype reference in FASTA format')
	# parser.add_argument('p_ref', help='P-type reference in FASTA format')
	parser.add_argument('ref', help='Reference sequence FASTA file')
	parser.add_argument('--genes', required=True, help='YAML of gene positions to be extracted')
	parser.add_argument('--gtype', action='store_true', help='YAML of gene positions to be extracted')
	parser.add_argument('--ptype', action='store_true', help='YAML of gene positions to be extracted')
	parser.add_argument('--prefix', required=True, help='YAML of gene positions to be extracted')
	parser.add_argument('--debug', action='store_true', help='YAML of gene positions to be extracted')

	return parser 

def init_aligner():
	aligner = Align.PairwiseAligner()
	aligner.mode = 'global'

	aligner.open_gap_score = -2.0
	aligner.extend_gap_score = -0.1
	aligner.target_end_gap_score = 0.0
	aligner.query_end_gap_score = 0.0
	
	return aligner

def extract_genes(aligner, ref, query, gene_dict, debug):

	# Perform pairwise alignment and get aligned sequences 
	align_iter = aligner.align(ref.seq, query.seq)
	ref_aligned, qry_aligned = next(align_iter)

	print(query.id)

	if debug:
		print(gene_dict)
		print("REFERENCE:")
		print(ref_aligned)
		print("QUERY:")
		print(qry_aligned)

	# Generate necessary alignment dictionaries to translate coordinates 
	align_dict, _ = make_align_dict(ref_aligned, qry_aligned)

	# STEP 1 - USE ALIGNMENT DICT TO TRANSLATE FROM REF POS TO ALIGNMENT POS
	# alignment_positions = align_dict[reference_positions]

	results = {}

	for gene, positions in gene_dict.items():
		align_start = align_dict[positions[0]]-1
		align_end = align_dict[positions[1]]-1

		if debug:
			print("GENE: " + gene)
			print("REFERENCE:")
			print(ref_aligned[align_start: align_end+1])
			print("QUERY:")
			print(qry_aligned[align_start: align_end+1])

		outseq = qry_aligned[align_start: align_end+1].replace("-","")

		outrecord = SeqIO.SeqRecord(
			Seq(outseq),
			id=query.id,
			name=query.name,
			description=query.description
		)

		results[gene] = outrecord

	return results

def main():
	parser = init_parser()
	args = parser.parse_args()

	# load reference gene positions from the args.genes YAML file 
	with open(args.genes,'r') as infile:
		gene_dict = yaml.safe_load(infile)

	aligner = init_aligner()

	# load the GI Norwalk reference 
	ref_seq = next(SeqIO.parse(args.ref, 'fasta'))

	# load the query sequence(s)
	qry_seq = list(SeqIO.parse(args.query, 'fasta'))

	if len(qry_seq) > 1:

		gene_results = list(map(lambda x : extract_genes(aligner, ref_seq, x, gene_dict, args.debug), qry_seq))
		#print(list(gene_results))
		gene_results = {gene: list(x[gene] for x in gene_results) for gene in gene_dict}

	else:
		gene_results = extract_genes(aligner, ref_seq, qry_seq[0], gene_dict, args.debug)
		#polymerase = extract_genes(aligner, qry_seq, ref_seq, gene_dict, debug)

	
	if args.gtype:
		SeqIO.write(gene_results['vp1'], args.prefix + '.vp1.fasta', 'fasta')
	if args.ptype:
		SeqIO.write(gene_results['rdrp'], args.prefix + '.rdrp.fasta', 'fasta')



if __name__ == '__main__':
	main()