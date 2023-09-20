import os 
import sys
import pysam
from collections import defaultdict
import re
import numpy as np
from Bio import Align


#%%
def parse_fasta(filepath):	
	seqs = {}
	with open(filepath, 'r') as handle:
		for line in handle.readlines():
			if line[0] == '>':
				header = line.strip().lstrip('>')
				seqs[header] = ''
			else:
				seqs[header] += line.strip()
	return seqs	

#%%
def write_fasta(seqs, outpath):
	try: 
		with open(outpath, 'w') as outfile:

			for header, seq in seqs.items():
				outfile.write(">" + header + '\n')
				outfile.write(seq + '\n')
	except Exception as e:
		print(str(e))
		return False
	return True

#%%
# helper function to handling read an input sam file
def alignment_reader(sam_input_path):
    """
    Function that returns a handle reading an input sam file
    """
    if sam_input_path:
        input_sam = pysam.AlignmentFile(sam_input_path, 'r')
    else:
        input_sam = pysam.AlignmentFile('-', 'r')
    return input_sam


def alignment_writer(bam_output_path, header):
    """
    Function that returns a handle that writes to a new output sam file
    """
    if bam_output_path:
        output_bam = pysam.AlignmentFile(bam_output_path, 'wb', header=header)
    else:
        output_bam = pysam.AlignmentFile('-', 'wb', header=header)

    return output_bam
#%%
def read_pair_generator(sam, region_string=None):
    """
    Generate read pairs for a SAM or BAM file or within a region string of said file.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in sam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


complement_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 
                    'W':'S', 'R':'Y', 'K':'M', 'Y':'R', 'S':'W', 'M':'K',
                    'B':'V', 'D':'H', 'H':'D', 'V':'B',
                    '*':'*', 'N':'N', '-':'-'}

def reverse_and_complement(seq):
    rseq = seq[::-1]
    rcseq = ''
    for i in rseq:  # reverse order
        rcseq += complement_dict[i]
    return rcseq

codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 
'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 
'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
'---':'-', 'XXX':'?'}

mixture_regex = re.compile('[WRKYSMBDHVN-]')

mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT', 'S':'CG', 
'M':'AC', 'V':'AGC', 'H':'ATC', 'D':'ATG', 
'B':'TGC', 'N':'ATGC', '-':'ATGC'}

ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.items())


def translate_nuc(seq, offset, resolve=False, return_list=False):
	"""
	Translate nucleotide sequence into amino acid sequence.
		offset by X shifts sequence to the right by X bases
	Synonymous nucleotide mixtures are resolved to the corresponding residue.
	Nonsynonymous nucleotide mixtures are encoded with '?' 
	"""
	
	seq = '-'*offset + seq
	
	aa_list = []
	aa_seq = ''	# use to align against reference, for resolving indels
	
	# loop over codon sites in nucleotide sequence
	for codon_site in range(0, len(seq), 3):
		codon = seq[codon_site:codon_site+3]
		
		if len(codon) < 3:
			break
		
		# note that we're willing to handle a single missing nucleotide as an ambiguity
		if codon.count('-') > 1 or '?' in codon:
			if codon == '---':	# don't bother to translate incomplete codons
				aa_seq += '-'
				aa_list.append(['-'])
			else:
				aa_seq += '?'
				aa_list.append(['?'])
			continue
		
		# look for nucleotide mixtures in codon, resolve to alternative codons if found
		num_mixtures = len(mixture_regex.findall(codon))
		
		if num_mixtures == 0:
			aa = codon_dict[codon]
			aa_seq += aa
			aa_list.append([aa])
			
		elif num_mixtures == 1:
			resolved_AAs = []
			for pos in range(3):
				if codon[pos] in mixture_dict.keys():
					for r in mixture_dict[codon[pos]]:
						rcodon = codon[0:pos] + r + codon[(pos+1):]
						if codon_dict[rcodon] not in resolved_AAs:
							resolved_AAs.append(codon_dict[rcodon])
							
			aa_list.append(resolved_AAs)
			
			if len(resolved_AAs) > 1:
				if resolve:
					# for purposes of aligning AA sequences
					# it is better to have one of the resolutions
					# than a completely ambiguous '?'
					aa_seq += resolved_AAs[0]
				else:
					aa_seq += '?'
			else:
				aa_seq += resolved_AAs[0]
				
		else:
			aa_seq += '?'
			aa_list.append(['?'])
			
	if return_list:
		return aa_list

	return aa_seq

def flex_translate(nt_seq, debug=False):
	
	min_frame = -1
	min_count = np.Inf
	best_seq = None

	for i in range(3):
		aa_seq = translate_nuc(nt_seq, i)
		aa_count = aa_seq.count("*")
		
		if debug: 
			print(aa_seq)
	
		if aa_count < min_count:
			min_count = aa_count
			best_seq = aa_seq
			min_frame = i
	
	# flip 1 and 2
	# this is to simplify the start position
	# instead of padding the start and shifting upstream (how this function does it), we just start further downstream instead, hence pad1 = ORF2 and vice versa
	if min_frame > 0:
		min_frame = 1 if min_frame == 2 else 2

	return best_seq, min_count, min_frame

def init_aligner(mode='global', open_gap=-1.0, x_gap=-0.1):
	aligner = Align.PairwiseAligner()
	aligner.mode = mode
	aligner.open_gap_score = open_gap
	aligner.extend_gap_score = x_gap
	aligner.target_end_gap_score = 0.0
	aligner.query_end_gap_score = 0.0
	return aligner


def make_align_dict(ref, qry):
	if len(ref) != len(qry):
		print("ERROR: Sequences are not aligned")
		return None
	align_dict = {} # {reference_position: alignment_position}
	qry_dict = {} # {alignment_position : qry_position}
	ref_idx = 0
	qry_idx = 0
	for idx, ref_char in enumerate(ref):
		if ref_char != "-":
			align_dict[ref_idx] = idx
			ref_idx += 1
		qry_dict[idx] = qry_idx
		if qry[idx] != '-':
			qry_idx += 1
	return align_dict, qry_dict

def get_boundaries (str):
    gap_prefix = re.compile('^[-]+')
    gap_suffix = re.compile('[-]+$')
    # return a tuple giving indices of subsequence without gap prefix and suffix
    res = [0,len(str)]
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        res[0] = len(left[0])
    if right:
        res[1] = len(str) - len(right[0])
    return res

if __name__ == '__main__':
	seqs = parse_fasta(sys.argv[1])

	for header, seq in seqs.items():
		print(">"+header)
		print(translate_nuc(seq, 0))
