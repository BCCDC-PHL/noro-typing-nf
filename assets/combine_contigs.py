from Bio import SeqIO
import os, sys
from glob import glob
import re 

contig_files = glob(os.path.join(sys.argv[1], "*.contigs.fa*"))

contigs = []

for f in contig_files:
	seqs = SeqIO.parse(f, 'fasta')
	name = re.split("[\._]", os.path.basename(f))[0]
	name = name.replace('-','')
	for s in seqs:
		s.id = name + "_" + s.id
		s.description = s.id
		s.name = s.id 

		contigs += (s,)

SeqIO.write(contigs, sys.argv[2], 'fasta')