# %%
import os, sys
import pandas as pd
import re
from Bio import SeqIO
from Bio import Entrez
from collections import defaultdict, Counter
import argparse

Entrez.email = os.environ['NCBI_EMAIL']
Entrez.api_key = os.environ['NCBI_API_KEY']

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('csv', help='Database headers in CSV format. Expected to contain < 5000 rows.')
	parser.add_argument('-a', '--accno_column', default='accno', help='Zero indexed position of the accession number column in the csv')
	parser.add_argument('-o', '--outfile', default="efetch.fasta", help='Scoring metric to select best reference/contig. Either bitscore (default), rawscore, or bsr.')
	return parser

def fetch_entries(accno_list, db_name, return_type, outpath):

	accnos = pd.Series(accno_list)

	accnos = accnos.drop_duplicates()
	
	print(accnos[0:10])

	query = Entrez.efetch(db=db_name, id=accnos, rettype=return_type, retmode='text')

	with open(outpath, 'w') as outfile:
		for i in query:
			outfile.write(i)
	
def parse_db_headers(filepath, accno_column):

	df = pd.read_csv(filepath)
	df = df.astype(str)
	df = df.apply(lambda x: x.str.replace(" ","_"))
	df = df.set_index(accno_column,drop=False)

	full_headers = df.apply(lambda x : "|".join(x), axis=1)

	return full_headers.to_dict()

def rename(seqs, rename_dict):

	for record in seqs:
		record.description = ''
		record.id = rename_dict[record.id.split(".")[0]]
	
	return seqs

def main():
	parser = get_parser()
	args = parser.parse_args()

	header_dict = parse_db_headers(args.csv, args.accno_column)

	# fetch_entries(accnos, 'nucleotide', 'fasta_cds_aa')
	fetch_entries(list(header_dict.keys()), 'nucleotide', 'fasta', args.outfile)

	seqs = list(SeqIO.parse(args.outfile, 'fasta'))

	#seqs = rename(seqs, header_dict)
	for s in seqs:
		s.description = s.id
		s.name = s.id 

	SeqIO.write(seqs, args.outfile, 'fasta')

if __name__ == '__main__':
	main()
