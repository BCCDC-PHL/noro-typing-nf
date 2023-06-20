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
	parser.add_argument('--header_pos_accno', default=0, help='Zero indexed position of accession number in FASTA header')
	parser.add_argument('--header_delim', default='|', help='Delimiter character used in the reference FASTA header')

	subparsers = parser.add_subparsers(help='sub-command help')
	# create the parser for the "a" command
	parser_a = subparsers.add_parser('database', help='One reference, multiple input queries. Extract genes from a BLAST DB of query inputs')
	parser_a.add_argument('-q', '--query', required=True, help='FASTA sequences to align to reference')
	parser_a.add_argument('-r', '--ref', required=True, help='Single reference sequence FASTA file. For norovirus, this is the G1 reference.')
	parser_a.add_argument('-p', '--positions', required=True, help='YAML file containing gene positions for extraction. Positions should be 1-indexed, as on GenBank, not 0-indexed.')
	parser_a.add_argument('-g', '--gene', required=True, help='Gene to extract.')
	parser_a.add_argument('-o', '--outfasta', help='Output FASTA file containing genes')
	parser_a.add_argument('-O', '--outyaml', required=True, help='Output YAML file with saved positions of all genes')
	parser_a.add_argument('-n', '--nthreads', default=1, type=int)
	parser_a.set_defaults(func=main_database)

	# create the parser for the "b" command
	parser_b = subparsers.add_parser('sample', help='Single input query mode. Extract one gene from one sample.')
	parser_b.add_argument('-q', '--query', required=True, help='Single FASTA sequence to align to chosen reference')
	parser_b.add_argument('-r', '--ref', required=True, help='Multi FASTA file of reference sequences')
	parser_b.add_argument('-p', '--positions', required=True, help='YAML file containing gene positions for extraction')
	parser_b.add_argument('-g', '--gene', required=True, help='Gene to extract. For norovirus, either vp1 or rdrp.')
	parser_b.add_argument('-o', '--outfasta', required=True, help='Output FASTA file containing genes')
	parser_b.set_defaults(func=main_sample)

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

def extract_gene_by_position(query, aligner, ref, gene_pos, debug):

	# Perform pairwise alignment and get aligned sequences 
	align_iter = aligner.align(ref.seq, query.seq)
	ref_aligned, qry_aligned = next(align_iter)

	align_start, align_stop = get_boundaries(ref_aligned)

	ref_aligned = ref_aligned[align_start: align_stop]
	qry_aligned = qry_aligned[align_start: align_stop]

	# Generate necessary alignment dictionaries to translate coordinates 
	align_dict, query_dict = make_align_dict(ref_aligned, qry_aligned)

	# STEP 1 - USE ALIGNMENT DICT TO TRANSLATE FROM REF POS TO ALIGNMENT POS
	align_start = align_dict[gene_pos[0]-1]
	align_end = align_dict[gene_pos[1]]  # JUSTIFIED BY SANITY CHECK 

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

	return outrecord, (query_dict[align_start], query_dict[align_end])

def sanity_check(seq):

	test_seq = copy(seq)

	if test_seq.seq.count("N") / len(test_seq) > 0.3:
		raise HighNContentException(f'ERROR: Query alignment has significant dropout in the target gene. Alignment failed. Exiting.\n{str(test_seq.seq)}')

	# remove excess bases that do not fit 
	if len(test_seq) % 3 != 0:
		extra_bases = len(test_seq) % 3
		test_seq.seq = test_seq.seq[:-extra_bases]


	aa_qry = test_seq.translate()

	if aa_qry.seq.count("*") > 1:
		print("WARNING: Translated nucleotide has multiple stop codons. Going to attempt flex translation.")
		print("REFERENCE: \n"+translate_nuc(str(seq.seq),0))
		print("QUERY: \n "+translate_nuc(str(seq.seq),0))
		
		seq_fix, stop_codon_count, frame_start = flex_translate(str(seq.seq))

		if stop_codon_count <= 1:
			print("FIXED.")
			print(seq_fix)
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

	# catch invalid gene inputs 
	if args.gene not in ['vp1', 'rdrp', 'all']:
		print("ERROR: Invalid gene selected.")
		sys.exit(1)
	
	# catch invalid filepath case 
	if args.gene == 'all' and (not re.search('\{[^\}]+\}', args.outfasta) or not re.search('\{[^\}]+\}', args.outyaml)):
		print("ERROR: Must include {gene} string in both the output FASTA/YAML paths when using gene = 'all' ")
		sys.exit(1)
	
	# parse gene inputs 
	if args.gene != 'all':
		genes = [args.gene]
	else:
		genes = ['vp1', 'rdrp']

	# iterate over all genes needing to be processed 
	for gene in genes:
		# instantiate a partial function used for multiprocessing
		mp_extract = partial(extract_gene_by_position, aligner=aligner, ref=ref_seq, gene_pos=position_dict[gene], debug=args.debug)

		# extract the gene from all query sequences
		with mp.Pool(args.nthreads) as pool:
			results = pool.map(mp_extract, qry_seqs)
		gene_seqs, gene_positions = zip(*results)
		
		# run a sanity check that fixes incorrect excisions
		gene_seqs = [sanity_check(x) for x in gene_seqs]

		# extract gene positions for later 
		gene_positions = { seq.id.split(args.header_delim)[args.header_pos_accno] : list([pos[0]+1, pos[1]+1]) for seq, pos in zip(qry_seqs, gene_positions)}

		# check for failed cases 
		failed = sum((1 for x in gene_seqs if not isinstance(x, SeqRecord)))

		# report failed cases 
		if failed > 0:
			print(f"WARNING: {failed} problematic sequences have been generated. \
				Only writing good quality ones. Sequence output will be lower than input.")
			
			gene_seqs = [x for x in gene_seqs if isinstance(x, SeqRecord)]

		# produce an output filename if none is given 
		if not args.outfasta:
			outpath = re.split("_|\.", os.path.basename(args.query))[0] + '_' + gene + '.fasta'
		elif args.gene == 'all':
			outpath = re.sub("\{[^\}]+\}", gene, args.outfasta)
			outyaml = re.sub("\{[^\}]+\}", gene, args.outyaml)
		else:
			outpath = args.outfasta
			outyaml = args.outyaml

		# write out the gene positions to a YAML file 
		with open(outyaml, 'w') as outfile:
			yaml.dump(gene_positions, outfile)
		
		# write output FASTA file 
		SeqIO.write(gene_seqs, outpath, 'fasta')

def main_sample(args):
	'''
	Main function for mode #1 (single query). Function used to extract a gene from a single query sequence 
	'''
	# catch invalid gene inputs 
	if args.gene not in ['vp1', 'rdrp', 'all']:
		print("ERROR: Invalid gene selected.")
		sys.exit(1)

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
	ref_seq_dict = {x.split(args.header_delim)[args.header_pos_accno]: y for x,y in ref_seq_dict.items()}

	if qry_seq.id.count("|") != 4:
		print("ERROR: Wrong sequence header input specified. Should contain 5 fields total: ID|TYPE|ACCNO|WORKFLOW|EXTRA")
		sys.exit(1)

	# parse out the reference accession number from the query FASTA file 
	accno_field = qry_seq.id.split("|")[2]

	if qry_seq.id.split("|")[4] == 'composite':

		accnos = accno_field.split("_")

		if (args.gene == 'vp1' and accnos[0] == 'NA') or (args.gene == 'rdrp' and accnos[1] == 'NA'):
			print("ERROR: Sequence failed in an earlier BLAST step. No reference accession found.")
			sys.exit(1)

		# this section can be expanded if more genes are desired
		if args.gene == 'vp1':
			ref_accno = accnos[0]
		elif args.gene == 'rdrp':
			ref_accno = accnos[1]
		else:
			print("ERROR: Not a valid gene entry.")
			sys.exit(1)

	else:
		ref_accno = accno_field

	# retrieve the appropriate reference using the chosen accession number
	ref_seq = ref_seq_dict[ref_accno]

	# load the positions from the YAML file 
	with open(args.positions,'r') as infile:
		position_dict = yaml.safe_load(infile)

	# extract the gene and perform sanity checks for the right reading frame 
	gene_record, _ = extract_gene_by_position(qry_seq, aligner, ref_seq, position_dict[ref_accno], args.debug)

	# catch failed cases where the N content is really high 
	try:
		gene_record = sanity_check(gene_record)
	except HighNContentException as e:
		print(e.message, file=sys.stderr)
		sys.exit(1)

	# check whether the gene is properly extracted 
	if gene_record:
		if not args.outfasta:
			output = re.split("_|\.", os.path.basename(args.query))[0]
		else:
			output = args.outfasta

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