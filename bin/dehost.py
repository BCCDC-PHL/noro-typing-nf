#!/usr/bin/env python3
import argparse
import csv
import pysam
import os
import sys
import subprocess
from tools import *

from collections import defaultdict
from pathlib import Path

def init_parser():
	'''
	Parser Arguments to pass to script from CL
	'''

	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--file', required=False, default=False, help='Composite Reference SAM file of mapped reads')
	parser.add_argument('-n', '--name', required=False, default=False, help='Name of the sample required when using stdin')
	parser.add_argument('-t', '--threads', required=False, type=int, default=1, help='Number of threads to run in parallel.')
	parser.add_argument('-r', '--ref_seqs', required=True, help='Reference ID of genome to keep. Default: MN908947.3')
	parser.add_argument('-q', '--keep_minimum_quality', required=False, type=int, default=60, help='Minimum quality of the reads to keep')
	parser.add_argument('-Q', '--remove_minimum_quality', required=False, type=int, default=0, help='Minimum quality of the reads to be included in removal')
	parser.add_argument('-o', '--output', required=False, default=False, help='Output BAM name')
	parser.add_argument('-R', '--revision', required=False, default='NA', help='Pass a pipeline commit hash to keep track of what version was ran')

	return parser


def remove_host_reads(samfile, keep_set, keep_minimum, remove_minimum):
	"""
	Function to remove human reads 
	Derived from the PHAC ncov-dehoster tool
	Filters reads based on the mapping qualities  
	"""
	header = samfile.header

	viral_reads = []
	human_count = 0
	unmapped_count = 0
	unmapped_reads = []
	poor_quality_count = 0

	for read_1, read_2 in read_pair_generator(samfile):
		# skip null reads 
		if read_1 is None or read_2 is None:
			continue

		# only look at primary alignments
		if (read_1.is_supplementary or read_1.is_secondary) and \
			(read_2.is_supplementary or read_2.is_secondary):
			continue

		if read_1.is_unmapped or read_2.is_unmapped:
			unmapped_count += 2
			unmapped_reads.append(read_1)
			unmapped_reads.append(read_2) 
			continue

		if (read_1.reference_name not in keep_set and read_1.mapping_quality >= remove_minimum) or \
			(read_2.reference_name not in keep_set and read_2.mapping_quality >= remove_minimum):
			human_count += 2
			
		elif read_1.mapping_quality >= keep_minimum and read_2.mapping_quality >= keep_minimum:
			viral_reads.append(read_1)
			viral_reads.append(read_2)

		else:
			poor_quality_count += 2

	return viral_reads, unmapped_reads, human_count, poor_quality_count

def main_phac():
	parser = init_parser()
	args = parser.parse_args()

	if args.file:
		sample_name = os.path.splitext(Path(args.file).stem)[0]
	elif not args.name:
		print("ERROR: Sample name not provided while using STDIN.")
		sys.exit(1)
	else:
		sample_name = args.name


	samfile = alignment_reader(args.file)

	read_list, human_filtered_count, poor_quality_count = remove_host_reads(samfile, args.keep_id, args.keep_minimum_quality, args.remove_minimum_quality )

	outfile = alignment_writer(args.output, template=samfile, header=samfile.header)

	for read in read_list:
		outfile.write(read)

	samfile.close()
	outfile.close()

	if len(read_list) == 0:
		percentage_kept = 0
	else:
		percentage_kept = len(read_list)/(human_filtered_count + poor_quality_count + len(read_list)) * 100

	line = {    'sample' : sample_name,
				'human_reads_filtered' : human_filtered_count, 
				'poor_quality_reads_filtered' : poor_quality_count,
				'paired_reads_kept' : len(read_list),
				'percentage_kept' : "{:.2f}".format(percentage_kept),
				'github_commit' : args.revision}
	
	with open(f'{sample_name}_stats.csv', 'w') as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames=line.keys())
		writer.writeheader()
		writer.writerow(line)



## ===============

def read_alignment(path):
	if not path:
		infile = pysam.AlignmentFile('-', 'r')
	elif path.endswith('.sam'):
		infile = pysam.AlignmentFile(path, 'r')
	else:
		infile = pysam.AlignmentFile(path, 'rb')
	return infile
	
def write_alignment(bam_output_path, header, mode='w'):
	"""
	Function that returns a handle that writes to a new output sam file
	"""
	if bam_output_path:
		output_bam = pysam.AlignmentFile(bam_output_path, mode, header=header)
	else:
		output_bam = pysam.AlignmentFile('-', mode, header=header)
	return output_bam

# Function to filter human reads
# Derived from ncov2019-artic-nf 
def filter_reads(input_sam_fp, output_bam_fp, contig_names):

	# use streams if args are None
	input_sam = read_alignment(input_sam_fp) 
	output_bam = write_alignment(output_bam_fp, header=input_sam.header)

	# if read isn't mapped or mapped to viral reference contig name
	viral_reads = 0
	human_reads = 0
	unmapped_reads = 0

	contig_dict = dict(zip(contig_names, [0] * len(contig_names)))

	# iterate over input from BWA
	for read in input_sam:
		# only look at primary alignments
		if not read.is_supplementary and not read.is_secondary:
			if read.reference_name in contig_dict:
				output_bam.write(read)
				contig_dict[read.reference_name] += 1
			elif read.is_unmapped:
				output_bam.write(read)
				unmapped_reads +=1
			else:
				human_reads += 1
	
	input_sam.close()
	output_bam.close()

	viral_reads = sum(contig_dict.values())
	total_reads = viral_reads + human_reads + unmapped_reads

	if total_reads > 0:
		print(f"viral read count = {viral_reads} "
			  f"({viral_reads/total_reads * 100:.2f}%)", file=sys.stderr)
		print(f"human read count = {human_reads} "
			  f"({human_reads/total_reads * 100:.2f}%) ", file=sys.stderr)
		print(f"unmapped read count = {unmapped_reads} "
			  f"({unmapped_reads/total_reads * 100:.2f}%)", file=sys.stderr)
		for contig in contig_dict:
			if contig_dict[contig] > 0:
				print(f'{contig} : {contig_dict[contig]}', file=sys.stderr)
	else:
		print(f"viral read count = {viral_reads} "
			  f"(0.00%)", file=sys.stderr)
		print(f"human read count = {human_reads} "
			  f"(0.00%) ", file=sys.stderr)
		print(f"unmapped read count = {unmapped_reads} "
			  f"(0.00%)", file=sys.stderr)

#%%
def main_artic():
	parser = init_parser()
	args = parser.parse_args()
	
	# if string is a file path, read the file and parse the set of reference sequences 
	if os.path.isfile(args.ref_seqs):
		print("Reading reference sequences from text file.")
		with open(args.ref_seqs, 'r') as infile:
			refs = [x.strip() for x in infile.readlines()]
	
	# if not a file, read string as the virus name 
	else:
		refs = list(str(args.ref_seqs))

	filter_reads(args.file, args.output, refs)

if __name__ == "__main__":
	main_artic()
