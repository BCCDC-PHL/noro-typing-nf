import os 
import sys
from Bio import SeqIO


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
def write_fasta(seq_dict, outpath):

	with open(outpath, 'w') as outfile:
		for header, seq in seq_dict.items():
			outfile.write('>' + header + '\n')
			outfile.write(seq + '\n')

def rename(inpath, fn, outpath):
	seqs = parse_fasta(inpath)

	for header in list(seqs.keys())[:]:
		newheader = fn(header)
		seqs[newheader] = seqs[header]
		del seqs[header]

	write_fasta(seqs, outpath)
	
	return "Complete"

def myfunction(header):
	if "(" in header:
		return header.replace("(",'-').replace(")","")
	else:
		return header
	
if __name__ == '__main__':
	res = rename(sys.argv[1], myfunction, sys.argv[2])
	