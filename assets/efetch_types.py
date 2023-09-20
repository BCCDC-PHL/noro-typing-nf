#!/usr/bin/env python3
# %%
import os, sys
import pandas as pd
from Bio import SeqIO, Entrez
from glob import glob
from datetime import datetime as dt
import argparse
import re

Entrez.email = os.environ['NCBI_EMAIL']
Entrez.api_key = os.environ['NCBI_API_KEY']

#%%
def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('fasta', help="Input FASTA file for which dates should be retrieved.")
	# parser.add_argument('-i', '--inplace', action='store_true', help="Indicates whether to append the dates to the end of the FASTA header with a | delimiter")
	parser.add_argument('-d', '--header-delim', default="|", type=str, help="Delimiter character used in the FASTA header.")
	parser.add_argument('-p', '--accno-pos', default=0, type=int, help="Zero-indexed position of the accession number in the FASTA header")
	parser.add_argument('-o', '--outfasta', help="Output FASTA with genotype/ptype fields appended.")
	parser.add_argument('-O', '--outtsv', help="Output genotype/ptype data in a metadata file. ")
	parser.add_argument('--debug', action='store_true')

	return parser.parse_args()

def num_to_roman(num):
	convert = {
		"1": "I",
		"2": "II",
		"3": "III",
		"4": "IV",
		"5": "V",
		"6": "VI",
		"7": "VII",
		"8": "VIII",
		"9": "IX",
		"10": "X",
	}
	return convert[str(num)]

def extract_type(string):
	"""
	Core worker function to resolve the wide variety of (messy) Genbank annotations
	Involved extensive testing on regex101.com
	"""
	main = re.compile('(.*[GP].*)[-_\/\[,\|](.*[GP].*)')
	no_ptype = re.compile("^[^P]+$", re.IGNORECASE)
	no_gtype = re.compile('^G[^PG\[]+P[^GP]+$', re.IGNORECASE)

	if string == '':
		return '', ''

	search1 = main.search(string)
	search2 = no_ptype.match(string)
	search3 = no_gtype.match(string)

	if search1:   # ideal case scenario; both genotype and ptype present 
		field1, field2 = main.search(string).groups()

		gtype, ptype = (field1, field2) if "P" in field2 else (field2, field1)

		# strip trailing whitespace and relevant characters
		gtype = gtype.strip()
		ptype = ptype.strip('] ')

		# fix the use of dash instead of period
		search = re.search('(G[IXV]+)-(.+)', gtype, re.IGNORECASE)
		if search:
			gtype = search.group(1) + '.' + search.group(2)

		# remove whitespace ; should not be whitespace by this point
		gtype = gtype.replace(" ", "-")
		ptype = ptype.replace(" ", "-")

		# g/ptype fixes to strain specific names
		if re.search("sydney", gtype, re.IGNORECASE):
			gtype = "GII.4-Sydney"
		if re.search("New.?Orleans", gtype, re.IGNORECASE): 
			gtype = "GII.4-New-Orleans"
		if re.search("\/\d+\/Changsha.+", gtype, re.IGNORECASE):
			gtype = re.sub("\/\d+\/Changsha.+", "", gtype)
		if re.search("Hong.?Kong", gtype, re.IGNORECASE):
			gtype = "GII.4-Hong-Kong"
		if re.search("New.?Orleans", ptype, re.IGNORECASE): 
			ptype = "GII.P4-New-Orleans"
		
		# remove extra bracketed annotation
		if re.search("\([^\)]+\)", ptype):
			ptype = re.sub("\([^\)]+\)","", ptype)

		# should not be multiple periods by this point
		# if multiple periods, delete them prior to the upcoming missing period fix 
		if ptype.count(".") > 1:
			print(gtype)
			print(ptype)
			ptype = ptype.replace(".","")
			
		# fix missing period in gtype 
		search = re.search("(G[IXV]+)[^\.\dP]?(\d+)", gtype)
		if search:
			gtype = search.group(1) + "." + search.group(2)

		# fix missing period in gtype, case 2
		search = re.search("(G[IXV]+)(N.+)", gtype)
		if search:
			gtype = search.group(1) + "." + search.group(2)

		# fix missing period in ptype 
		search = re.search("(G[IXV]+)(P.+)", ptype)
		if search:
			if "." in ptype:
				ptype = ptype.replace(".","")
			ptype = search.group(1) + "." + search.group(2)
	
		# fix for an oversimplified ptype call 
		if re.search("^P.+", ptype):
			gbase = re.search("(G[IXV]+)\.", gtype)

			if gbase: 
				ptype = gbase.group(1) + "." + ptype
			else:
				print("Failed to fix abbrev P type ")
		
		# fix for extra period within the ptype call
		search = re.search('(P)[IXV]*\.(\d+)', ptype)
		if search:
			ptype = search.group(1) + search.group(2)

		return gtype, ptype
	

	elif search2:     # secondary case with only genotype present 
		if not string.startswith("G"):
			string = "G"+string

		# handle case of numeric value 
		search = re.search("G(\d+)", string)
		if search:
			string = "G" + num_to_roman(search.group(1))

		string = string.replace("-",'.')
		string = string.replace("_","-")
		return string, ""

	elif search3:  	 # tertiary case with only Ptype present 
		string = string.replace("_","-")

		# fix missing period in ptype 
		search = re.search("(G[IXV]+)(P.+)", string)
		if search:
			string = search.group(1) + "." + search.group(2)

		return "", string

	else:			# failed to resolve to one of the three categories above
		return "FAIL", string
	
		
def fetch_types(accno_list):

	if len(accno_list) == 0:
		print("ERROR: At least one accession number is required to query Genbank")
		sys.exit(1)
		
	query = Entrez.efetch(db='nucleotide', id=accno_list, rettype='gb', retmode='text')

	seq_types = {}

	expr_accno = re.compile('ACCESSION\s+([^\s]+)')
	expr_type = re.compile('genotype:\s(.+)\"')

	accno = ''
	for line in query:
		line = line.strip("\n ")

		if line.startswith("ACCESSION"):
			# start of the next gb file 
			print(line)
			accno = expr_accno.search(line).group(1)
			seq_types[accno] = ''

		# retrieve collection dates (preferred date choice)
		if line.startswith('/note="genotype'):
			if expr_type.search(line): 
				stype = expr_type.search(line).group(1)
				seq_types[accno] = stype
			else:
				print(line)

	return seq_types

#%%
def main():
	args = init_parser()

	# parse sequences in FASTA format 
	seqs = list(SeqIO.parse(args.fasta, 'fasta'))

	# retrieve the accession numbers from the FASTA file 
	for seq in seqs:
		accno = seq.id.split(args.header_delim)[args.accno_pos].split(".")[0]
		seq.id = accno
		seq.name = accno
		seq.description = accno

	# retrieve the sequence types  
	if os.path.isfile(".types.cache"):
		with open(".types.cache", 'r') as infile:
			fields = [x.strip().split(",") for x in infile.readlines()]
			stype_dict = {x[0] : x[1] for x in fields}
	else:
		stype_dict = fetch_types([s.id for s in seqs])
		with open(".types.cache", 'w') as outfile:
			for name, stype in stype_dict.items(): 
				outfile.write(name + "," + stype + "\n")

	# fix sequence types into a common format 
	stype_fixed = {name.split(".")[0]: extract_type(stype) for name, stype in stype_dict.items()}

	print("Excluding empty typing sequences")
	print(f"Before: {len(seqs)} {len(stype_fixed)}")
	seqs = [x for x in seqs if stype_fixed[x.id][0] != '' or stype_fixed[x.id][1] != '' ]
	stype_fixed = {x:y for x, y in stype_fixed.items() if y[0] != '' or y[1] != ''}
	print(f"After: {len(seqs)} {len(stype_fixed)}")
	
	# rename the FASTA sequences 
	for seq in seqs:
		newheader = "|".join([seq.id, "_".join(stype_fixed[seq.id]), 'global'])
		seq.id = newheader
		seq.name = newheader
		seq.description = newheader

	# reformat into a dataframe 
	stype_rows = [(x, y[0], y[1]) for x, y in stype_fixed.items()]
	stype_df = pd.DataFrame(stype_rows, columns=['accno','genotype','ptype'])

	# output TSV metadata 
	stype_df.to_csv(args.outtsv, sep='\t', index=False)

	# output FASTA file 
	SeqIO.write(seqs, args.outfasta, 'fasta')

#%%
if __name__ == '__main__':
	main()
# %%
