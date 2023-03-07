#!/usr/bin/env python3
import argparse
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument('reference', help='Reference file in FASTA format')
parser.add_argument('bam', help='BAM file of reads aligned to reference' )

args = parser.parse_args()
bamfile = pysam.AlignmentFile(args.bam, "rb")
reference = pysam.FastaFile(args.reference)

print("\t".join(["contig", "position", "depth", "ref_base", "count_A", "count_C", "count_G", "count_T", "count_del", "alt_frequency"]))
for plup_column in bamfile.pileup(ignore_orphans=False):
    freqs = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0, 'R': 0 }
    for plup_read in plup_column.pileups:
        if plup_read.is_del:
            freqs['-'] += 1
        elif plup_read.is_refskip:
            freqs['R'] += 1
        else:
            base = plup_read.alignment.query_sequence[plup_read.query_position]
            freqs[base] += 1

    reference_base = reference.fetch(plup_column.reference_name, plup_column.reference_pos, plup_column.reference_pos+1)

    # count alt depth
    bc_depth = 0
    alt_depth = 0
    for base in "ACGT":
        bc_depth += freqs[base]
        if base != reference_base:
            alt_depth += freqs[base] 

    if bc_depth > 0:
        alt_freq = alt_depth / bc_depth
    else:
        alt_freq = 0.0

    out = list()

    outputs = [
        plup_column.reference_name,
        plup_column.pos + 1, # output 1-based coordinates for consistency with samtools mpileup
        bc_depth + freqs['-'],
        reference_base,
        freqs['A'],
        freqs['C'],
        freqs['G'],
        freqs['T'],
        freqs['-'],
        alt_freq
    ]
    outputs = [str(x) for x in outputs]
    
    print("\t".join(outputs))
