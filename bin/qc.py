#!/usr/bin/env python3

from Bio import SeqIO
import csv
import subprocess
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import shlex
import numpy as np

"""
This script can incorporate as many QC checks as required
as long as it outputs a csv file containing a final column
headed with 'qc_pass' and rows for each sample indcating
'TRUE' if the overall QC check has passed or 'FALSE' if not.
"""

def make_qc_plot(depth_pos, n_density, samplename, window=200):
    depth_df = pd.DataFrame( { 'position' : [pos[1] for pos in depth_pos], 'depth' : [dep[2] for dep in depth_pos] } )
    depth_df['depth_moving_average'] = depth_df.iloc[:,1].rolling(window=window).mean()

    n_df = pd.DataFrame( { 'position' : [pos[0] for pos in n_density], 'n_density' : [dens[1] for dens in n_density] } )

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel('Position')

    ax1.set_ylabel('Depth', color = 'g')
    ax1.set_ylim(top=10**5, bottom=1)
    ax1.set_yscale('log')
    ax1.plot(depth_df['depth_moving_average'], color = 'g')

    ax2.set_ylabel('N density', color = 'r')  
    ax2.plot(n_df['n_density'], color = 'r')
    ax2.set_ylim(top=1)

    plt.title(samplename)
    plt.savefig(samplename + '.depth.png')

def read_depth_file(bamfile):
    p = subprocess.Popen(['samtools', 'depth', '-a', bamfile],
                       stdout=subprocess.PIPE)
    out, err = p.communicate()
    counter = 0

    pos_depth = []
    for ln in out.decode('utf-8').split("\n"):
       if ln:
          pos_depth.append(ln.split("\t"))
    
    return pos_depth


def get_covered_pos(pos_depth, min_depth):
    counter = 0
    for contig, pos,depth in pos_depth:
        if int(depth) >= min_depth:
            counter = counter + 1
    
    return counter

def get_N_positions(fasta):
    n_pos =  [i for i, letter in enumerate(fasta.seq.lower()) if letter == 'n']

    return n_pos

def get_pct_N_bases(fasta):
    
    count_N = len(get_N_positions(fasta))

    pct_N_bases = count_N / len(fasta.seq) * 100

    return pct_N_bases

def get_largest_N_gap(fasta):
    n_pos = get_N_positions(fasta)

    n_pos = [0] + n_pos + [len(fasta.seq)]

    n_gaps = [j-i for i, j in zip(n_pos[:-1], n_pos[1:])]

    return sorted(n_gaps)[-1]

def get_ref_length(ref):
    record = SeqIO.read(ref, "fasta")
    return len(record.seq)

def get_depth_quantiles(pos_depth):
    _ , _ , depths = zip(*pos_depth)
    depths = [int(x) for x in depths]
    depths = np.array(depths)

    quantiles = [
        np.min(depths),
        *np.quantile(depths, [0.1, 0.25, 0.50]).tolist(),
        np.max(depths),
        np.mean(depths)
    ]

    return quantiles

def sliding_window_N_density(sequence, window=10):

    sliding_window_n_density = []
    for i in range(0, len(sequence.seq), 1):
        window_mid = i + ( window / 2)
        window_seq = sequence.seq[i:i+window]
        n_count = window_seq.lower().count('n')
        n_density = n_count / window

        sliding_window_n_density.append( [ window_mid, n_density ] )

    return sliding_window_n_density

def get_num_reads(bamfile):

    st_filter = '0x900'
    command = 'samtools view -c -F{} {}'.format(st_filter, bamfile)
    what = shlex.split(command)

    return subprocess.check_output(what).decode().strip()
    
def run_qc(args):
    ## Depth calc
    ref_length = get_ref_length(args.ref)
    depth_pos = read_depth_file(args.bam)

    depth_covered_bases = get_covered_pos(depth_pos, args.min_depth)

    pct_covered_bases = depth_covered_bases / ref_length * 100

    ## Number of aligned reads calculaton
    num_reads = get_num_reads(args.bam)

    # Unknown base calcs
    fasta = SeqIO.read(args.consensus, "fasta")

    pct_N_bases   = 0
    largest_N_gap = 0
    qc_pass       = "FALSE"

    if len(fasta.seq) != 0:

        pct_N_bases = get_pct_N_bases(fasta)
        largest_N_gap = get_largest_N_gap(fasta)

    	# QC PASS / FAIL
        if pct_N_bases < 50.0 and pct_covered_bases > 0.80:
            qc_pass = "TRUE"

    quantiles = get_depth_quantiles(depth_pos)

    qc_line = { 'sample_name' : args.sample_name,
                'depth_min'    : "{:.2f}".format(quantiles[0]),
                'depth_10'    : "{:.2f}".format(quantiles[1]),
                'depth_25'    : "{:.2f}".format(quantiles[2]),
                'depth_50'    : "{:.2f}".format(quantiles[3]),
                'depth_max'    : "{:.2f}".format(quantiles[4]),
                'depth_mean'   : "{:.2f}".format(quantiles[5]),
                'pct_N_bases' : "{:.2f}".format(pct_N_bases),
                f'pct_covered_{args.min_depth}' : "{:.2f}".format(pct_covered_bases), 
                'longest_no_N_run' : largest_N_gap,
                'num_aligned_reads' : num_reads,
                #    'fasta': args.consensus, 
                #     'bam' : args.bam,
                'qc_pass' : qc_pass}

    df = pd.DataFrame(qc_line, index=[0])
    df.to_csv(args.outfile, index=False)

    N_density = sliding_window_N_density(fasta)
    make_qc_plot(depth_pos, N_density, args.sample_name)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', required=True)
    parser.add_argument('-s', '--sample_name', required=True)
    parser.add_argument('-r', '--ref', required=True)
    parser.add_argument('-b', '--bam', required=True)
    parser.add_argument('-c', '--consensus', required=True)
    parser.add_argument('-d', '--min-depth', default=10, type=int)

    args = parser.parse_args()
    run_qc(args)

if __name__ == "__main__":
    main()
