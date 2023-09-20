from Bio import SeqIO, Align, AlignIO
import sys, os
import argparse
from collections import Counter, defaultdict

# usage: "python align2alleles.py --reference-name MN908947.3 "qc_analysis/{prefix}_aligned.fasta" > {output}"
# produces an alleles.tsv output file
# write results in matrix form where rows are samples
# and columns are variant positions
def write_result_matrix(variant_positions, sequences):
    # print header
    # +1 is to report 1-based coordinates
    print("\t".join(["strain"] + [str(i + 1) for i in variant_positions]))

    for s in sequences:
        out = list()
        out.append(s[0])

        for i in variant_positions:
            out.append(str(s[1][i]))

        print("\t".join(out))

# write results as a TSV file with one row per variant found in a sample

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--min-allele-count', default=1, type=int)
    parser.add_argument('--mode', default="variant_list", type=str)
    parser.add_argument('-r', '--reference', required=True, help="Reference sequence in FASTA format")
    parser.add_argument('-q', '--query', required=True, help="Query sequence in FASTA format")
    return parser.parse_args()


def main():
    args = parse_arguments()

    ref_fasta = next(SeqIO.parse(args.reference, 'fasta'))
    query_fasta = next(SeqIO.parse(args.query, 'fasta'))

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'

    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    alignment = next(aligner.align(ref_fasta.seq, query_fasta.seq))

    align_ref, align_query = alignment

    align_chars = zip(align_ref, align_query)

    # # discover variant columns
    # counters = list()
    # for i in range(0, alen):
    #     counters.append(collections.Counter())

    # #
    # reference_sample = None

    # # count the number of occurrences of each base at each position
    # for s in sequences:

    #     if s.name == args.reference_name:
    #         reference_sample = s

    #     for i in range(0, alen):
    #         b = s.sequence[i]
    #         counters[i].update(b)

    # if reference_sample is None:
    #     sys.stderr.write("error: reference sample could not be found")
    #     sys.exit(1)

    # write results

    # determine which positions have a SNP
    variants = []
    for n, (ref_char, qry_char) in enumerate(align_chars, start=1):
        if ref_char != 'N' and qry_char != 'N' and ref_char != qry_char:
            variants.append(n)


    # def write_variant_list(variant_positions, counters, sequences, reference_sample):
    #     print("\t".join(["pos", "ref_allele", "alt_allele"]))

    #     for s in sequences:
    #         for i in variant_positions:
    #             b = s.sequence[i]
    #             if b != reference_sample.sequence[i]:
    #                 print("\t".join([str(i+1), reference_sample.sequence[i], b]))

    ambiguous_counts = defaultdict(int)
    ambiguous_alleles = defaultdict(dict)
    for sample, records in alleles.data.items():
        for position in records:
            aa = records[position]['alt']
            if aa != 'N' and ncov.parser.is_variant_iupac(aa):
                p = int(position)
                ambiguous_counts[p] += 1
                ambiguous_alleles[p][aa] = 1

    print("position\tcount\talleles")
    for position in sorted(ambiguous_counts.keys()):
        count = ambiguous_counts[position]
        if count >= args.min_count:
            print("%d\t%d\t%s" % (int(position), count, ",".join(ambiguous_alleles[position])))

    # elif args.mode == "variant_frequency":
    #     for i in range(0, alen):
    #         reference_base = reference_sample.sequence[i]
    #         for b in ("A", "C", "G", "T"):
    #             if b == reference_base:
    #                 continue
    #             c = counters[i][b]
    #             if c > 0:
    #                 out = [ reference_name, str(i + 1), reference_base, b, counters[i][reference_base], counters[i][b] ]
    #                 print("\t".join([str(x) for x in out]))

if __name__ == '__main__':
    main()

