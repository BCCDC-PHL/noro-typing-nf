#!/usr/bin/env python3
#%%
import os, sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
import re
import shutil

#%%
def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('bam_file', help='Input BAM file')
	parser.add_argument('ref_fasta',  help='Final pipeline references in FASTA format')
	#parser.add_argument('pfasta',  help='P-type reference in FASTA format')
	parser.add_argument('--outbam', required=True, help='Merged BAM output file.')
	parser.add_argument('--outfasta', required=True, help='FASTA of the new merged reference sequence.')
	# NOTE: the following arguments should never change assuming "filter_fasta.py" has been run
	parser.add_argument('--header_pos_accno', default=0, help='Zero indexed position of accession number')
	parser.add_argument('--header_pos_type', default=1, help='Zero indexed position of sequence type') 
	parser.add_argument('--header_delim', default='|', help='Delimiter for incoming sequence headers')  
	return parser.parse_args()


def make_merged_fasta(bampath, fastapath, combined_name, debug=False):

	infile = pysam.AlignmentFile(bampath, 'rb')

	ref_dict = SeqIO.to_dict(SeqIO.parse(fastapath, 'fasta'))

	if len(ref_dict) != 2:
		print("ERROR: Single reference sequence not provided for either genotype or ptype FASTA.")

	# remove later
	# ref_dict.update((x.replace("-","_"), y) for x, y in list(ref_dict.items()))

	combined_ref = {}
	seq_len = 0

	for pl_column in infile.pileup():
		pos = pl_column.reference_pos
		ref_name = pl_column.reference_name
		seq_len = max(seq_len, pos)

		if pos in combined_ref:
			if debug: 
				print(pos)
				print(ref_name)
				print(pl_column.get_num_aligned())
			
			if combined_ref[pos][1] < pl_column.get_num_aligned():
				combined_ref[pos] = (str(ref_dict[ref_name].seq[pos: pos+1]), pl_column.get_num_aligned(), pl_column.reference_name)
		else:
			combined_ref[pos] = (str(ref_dict[ref_name].seq[pos: pos+1]), pl_column.get_num_aligned(), pl_column.reference_name)
	
	seq_len += 1  # position values are 0-indexed; this corrects seq_len to be 1-indexed

	print("Filling in empty positions.")
	
	empty_pos = 0 
	for i in range(seq_len):
		if i not in combined_ref:
			empty_pos += 1
			combined_ref[i] = ('N', 0, 'NA')
	
	print(f"Empty positions: {empty_pos}")

	bases = [(x, y[0]) for x, y in combined_ref.items()]
	bases = sorted(bases, key=lambda x: x[0])
	sequence = "".join(x[1] for x in bases)

	return SeqRecord(
		Seq(sequence),
		id=combined_name,
		name=combined_name,
		description=combined_name
	)


def make_merged_bam(bampath, new_ref_name, new_ref_len, outbam):

	bamfile = pysam.AlignmentFile(bampath, 'rb')

	# modify the original BAM header to include the new combined name
	new_header = bamfile.header.to_dict()
	new_header['SQ'].append({'SN': new_ref_name, 'LN': new_ref_len})

	print("Writing to temporary BAM...")
	with pysam.AlignmentFile('temp.sam', 'w', header=new_header) as outfile:
		for n, read in enumerate(bamfile.fetch()):
			# if n == 0:
			# 	temp_name = read.reference_name
			# read.reference_name = temp_name
			outfile.write(read)
	bamfile.close()
	
	# print("Sorting temp.bam...")
	# pysam.sort("-o", "temp.sorted.bam", "temp.bam")
	# print("Indexing temp.sorted.bam...")
	# pysam.index('temp.sorted.bam')

	print("Writing final BAM output with new header & read names...")
	with pysam.AlignmentFile('temp.sam', 'r') as infile, \
		pysam.AlignmentFile('temp.bam', 'wb', header=new_header) as outfile:
		for read in infile.fetch():
			read.reference_name = new_ref_name
			outfile.write(read)
	
	print("Sorting temp.bam...")
	pysam.sort("-o", outbam, "temp.bam")
	print("Indexing final output BAM file...")
	pysam.index(outbam)

	print("Cleaning up temp files...")
	os.remove('temp.sam')
	os.remove('temp.bam')
	# os.remove('temp.sorted.bam')
	# os.remove('temp.sorted.bam.bai')

	print("Complete.")


#%%
def main():
	args = get_args()

	delim = args.header_delim
	a_pos = args.header_pos_accno
	t_pos = args.header_pos_type

	bamfile = pysam.AlignmentFile(args.bam_file, 'rb')
	N_ref = bamfile.nreferences
	bamfile.close()

	if N_ref > 1:
		fasta = list(SeqIO.parse(args.ref_fasta,'fasta'))

		# create a combined reference name 
		# ref_names = [bamfile.get_reference_name(i).split(delim) for i in range(bamfile.nreferences)]
		# accnos = [x[t_pos].replace("_","-") for x in ref_names]
		# types = [x[a_pos].replace("_",'-') for x in ref_names]
		# combined_ref_name = "_".join(accnos) + "|" + "_".join(types)

		if len(fasta[0].description.split()) < 2:
			print("ERROR: Missing full header as the description entry.")
			sys.exit(1)

		combined_ref_name = fasta[0].description.split()[-1]

		# compute the maximum sequence length 
		combined_ref_length = max((len(seq) for seq in fasta)) 

		outfasta = make_merged_fasta(args.bam_file, args.ref_fasta, combined_ref_name)

	elif N_ref == 1:
		outfasta = next(SeqIO.parse(args.ref_fasta,'fasta'))

		if len(outfasta.description.split()) < 2:
			print("ERROR: Missing full header as the description entry.")
			sys.exit(1)
		
		# ref_fields = bamfile.get_reference_name(0).split(delim)
		# combined_ref_name = "_".join([ref_fields[a_pos]]*2) + "|" + "_".join(ref_fields[t_pos]*2)
		combined_ref_name = outfasta.description.split()[-1]
		combined_ref_length = len(outfasta)

		outfasta.id = combined_ref_name
		outfasta.name = combined_ref_name
		outfasta.description = combined_ref_name

		# shutil.copy2(args.bam_file, args.outbam)
		# pysam.index(args.outbam)
	else:
		print("ERROR: Looks like a problematic BAM file input.")
		sys.exit(1)

	make_merged_bam(args.bam_file, combined_ref_name, combined_ref_length, args.outbam)

	# if there are two different references, 
	# create a combined reference FASTA by taking the winning base at every position 
	SeqIO.write(outfasta, args.outfasta, 'fasta')
	print("Finished writing FASTA output.")

#%%
if __name__ == '__main__':
	main()
