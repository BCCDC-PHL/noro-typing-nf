# %%
import os, sys
import pandas as pd
import re
from Bio import SeqIO
from Bio import Entrez
from collections import defaultdict, Counter

Entrez.email = os.environ['NCBI_EMAIL']
Entrez.api_key = os.environ['NCBI_API_KEY']

def fetch_entries(accno_list, db_name, return_type):

	accnos = pd.Series(accno_list)

	accnos = accnos.drop_duplicates()
	
	print(accnos.to_list()[0:10])

	query = Entrez.efetch(db=db_name, id=accnos.to_list(), rettype=return_type, retmode='text')

	with open('efetch.fasta', 'w') as outfile:
		for i in query:
			outfile.write(i)


if __name__ == '__main__':

	if len(sys.argv) > 1:
		with open(sys.argv[1], 'r') as handle:
			accnos = [x.lstrip(">").strip() for x in handle.readlines()]
			print(accnos[0:10])
	else:
		# simple test case
		accnos = ['BBB94175.1']

	fetch_entries(accnos, 'nucleotide', 'fasta_cds_aa')


