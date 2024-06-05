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
from functools import partial
import multiprocessing as mp

def init_parser():
	# create the top-level parser
	parser = argparse.ArgumentParser(prog='PROG')
	parser.add_argument('--debug', action='store_true')
	parser.add_argument('--ref_header_pos_accno', default=0, help='Zero indexed position of accession number in FASTA header. Default 0 (1st)')
	parser.add_argument('--ref_header_delim', default='|', help='Delimiter character used in the reference FASTA header')
	parser.add_argument('-q', '--query', required=True, help='MultiFASTA of sequences from which genes will be extracted.')
	parser.add_argument('-r', '--ref', required=True, help='Single reference sequence FASTA file. For norovirus, this is the G1 reference.')
	parser.add_argument('-p', '--positions', required=True, help='YAML file containing gene positions for extraction. Positions should be 1-indexed, as on GenBank, not 0-indexed.')
	parser.add_argument('-g', '--gene', required=True, help='Name of the gene to extract.')
	parser.add_argument('-o', '--outfasta', required=True, help='Output FASTA file containing genes')
	parser.add_argument('-O', '--outyaml', required=True, help='Output YAML file with saved positions of all genes')
	parser.add_argument('-n', '--nthreads', default=1, type=int)

	return parser 

class HighNContentException(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def init_aligner():
	aligner = Align.PairwiseAligner()
	aligner.mode = 'global'
	aligner.open_gap_score = -2.0
	aligner.extend_gap_score = -0.1
	aligner.target_end_gap_score = 0.0
	aligner.query_end_gap_score = 0.0
	
	return aligner

def clean_sequence(seq):

	record = copy(seq)

	if record.seq.count("N") / len(record) > 0.3:
		raise HighNContentException(f'ERROR: Query alignment has significant dropout in the target gene. Alignment failed. Exiting.\n{str(record.seq)}')

	# remove excess bases that do not fit 
	if len(record) % 3 != 0:
		extra_bases = len(record) % 3
		record.seq = record.seq[:-extra_bases]


	aa_qry = record.translate()

	if aa_qry.seq.count("*") > 1:
		print("WARNING: Translated nucleotide has multiple stop codons. Attempting flex translation.")
		print("REFERENCE: \n"+translate_nuc(str(record.seq),0))
		print("QUERY: \n "+translate_nuc(str(record.seq),0))
		
		seq_fix, stop_codon_count, frame_start = flex_translate(str(record.seq))

		if stop_codon_count <= 1:
			print("FIXED.")
			print(seq_fix)
			record.seq = record.seq[frame_start:]
			return record
			
		else:
			print("WARNING: Problematic sequence.")
			print("REFERENCE: \n"+translate_nuc(str(record.seq),0))
			print("QUERY: \n"+translate_nuc(str(record.seq),0))

			return None

	return record

def extract_gene(query, aligner, ref, gene_pos, debug):

	# Perform pairwise alignment and get aligned sequences 
	ref_aligned, qry_aligned = next(aligner.align(ref.seq, query.seq))

	align_start, align_stop = get_boundaries(ref_aligned)

	ref_aligned = ref_aligned[align_start: align_stop]
	qry_aligned = qry_aligned[align_start: align_stop]

	# Generate necessary alignment dictionaries to translate coordinates 
	align_dict, query_dict = make_align_dict(ref_aligned, qry_aligned)

	# STEP 1 - USE ALIGNMENT DICT TO TRANSLATE FROM REF POS TO ALIGNMENT POS
	align_start = align_dict[gene_pos['start']-1]
	align_end = align_dict[gene_pos['end']]  # JUSTIFIED BY SANITY CHECK 

	print(query.id)

	# formulate the final gene extracted from the query 
	outseq = qry_aligned[align_start: align_end].replace("-","")

	if debug:
		print(f"Reference Length: {len(ref)}")
		print(f"Query Length: {len(query)}")
		print(f"Output Length: {len(outseq)} % 3 == {len(outseq) % 3}")
		print("REFERENCE:")
		print(ref_aligned)
		print("QUERY:")
		print(qry_aligned)
		print("REFERENCE:")
		print(ref_aligned[align_start: align_end+1])
		print("QUERY:")
		print(qry_aligned[align_start: align_end+1])
		print(f"Alignment Position Start: {align_start}")
		print(f"Alignment Position End: {align_end}")
		print(f"Query Position Start: {query_dict[align_start]}")
		print(f"Query Position End: {query_dict[align_end]}")
		print("TEST")
		print(query.seq[query_dict[align_start]:query_dict[align_end]])
		print(outseq)
		print("REFERENCE:")
		print(translate_nuc(ref_aligned[align_start: align_end+1],0))
		print("QUERY:")
		print(translate_nuc(qry_aligned[align_start: align_end+1].replace("-",""),0))
		print(f"Output Length: {len(outseq)} % 3 == {len(outseq) % 3}")

	outrecord = SeqIO.SeqRecord(
		Seq(outseq),
		id=query.id,
		name=query.name,
		description=query.description
	)

	outrecord = clean_sequence(outrecord)

	return outrecord, (query_dict[align_start], query_dict[align_end])


def main(args):
	print("Extracting genes from database...")

	# load reference gene positions from the args.genes YAML file 
	with open(args.positions,'r') as infile:
		position_dict = yaml.safe_load(infile)

	aligner = init_aligner()

	# load the reference seq list 
	ref_seq = next(SeqIO.parse(args.ref, 'fasta'))

	# load the database (a collection of query sequences)
	qry_seqs = SeqIO.parse(args.query, 'fasta')

	# catch invalid gene inputs 
	if args.gene not in list(position_dict.keys()) + ['all']:
		print("ERROR: Invalid gene selected.")
		sys.exit(1)
	
	# catch invalid filepath case 
	if args.gene == 'all' and (not re.search(r'\{[^\}]+\}', args.outfasta) or not re.search(r'\{[^\}]+\}', args.outyaml)):
		print("ERROR: Must include {gene} string in both the output FASTA/YAML paths when using gene = 'all' ")
		sys.exit(1)
	
	# parse gene inputs 
	if args.gene != 'all':
		genes = [args.gene]
	else:
		genes = list(position_dict.keys())

	# iterate over all genes needing to be processed 
	for gene in genes:
		# instantiate a partial function used for multiprocessing
		mp_extract_fn = partial(extract_gene, aligner=aligner, ref=ref_seq, gene_pos=position_dict[gene], debug=args.debug)

		# extract the gene from all query sequences
		with mp.Pool(args.nthreads) as pool:
			results = pool.map(mp_extract_fn, qry_seqs)
		gene_seqs, gene_positions = zip(*results)

		# check for failed cases 
		failed = sum((1 for x in gene_seqs if not isinstance(x, SeqRecord)))

		# report failed cases 
		if failed > 0:
			print(f"WARNING: {failed} problematic sequences have been generated. \
				Only writing good quality ones. Sequence output will be lower than input.")
			
			gene_seqs = [x for x in gene_seqs if isinstance(x, SeqRecord)]

		# extract gene positions for later 
		gene_positions = { seq.id.split(args.ref_header_delim)[args.ref_header_pos_accno] : [pos[0]+1, pos[1]+1] for seq, pos in zip(qry_seqs, gene_positions)}

		# produce an output filename if none is given 
		if not args.outfasta:
			outpath = re.split(r"_|\.", os.path.basename(args.query))[0] + '_' + gene + '.fasta'
		elif args.gene == 'all':
			outpath = re.sub(r"\{[^\}]+\}", gene, args.outfasta)
			outyaml = re.sub(r"\{[^\}]+\}", gene, args.outyaml)
		else:
			outpath = args.outfasta
			outyaml = args.outyaml

		# write out the gene positions to a YAML file 
		with open(outyaml, 'w') as outfile:
			yaml.dump(gene_positions, outfile)
		
		# write output FASTA file 
		SeqIO.write(gene_seqs, outpath, 'fasta')

if __name__ == '__main__':
	args = init_parser().parse_args()

	# call the assigned default function chosen by the user input 
	main(args)