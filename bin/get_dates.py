#!/usr/bin/env python3
# %%
import os, sys
import pandas as pd
from Bio import SeqIO, Entrez
from glob import glob
from datetime import datetime as dt
import argparse
import re

#%%
def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('fasta', help="Input FASTA file for which dates should be retrieved.")
	# parser.add_argument('-i', '--inplace', action='store_true', help="Indicates whether to append the dates to the end of the FASTA header with a | delimiter")
	parser.add_argument('-d', '--header-delim', default="|", type=str, help="Delimiter character used in the FASTA header.")
	parser.add_argument('-p', '--accno-pos', default=0, type=int, help="Zero-indexed position of the accession number in the FASTA header")
	parser.add_argument('-o', '--output', help="Output dates file name. Filename chosen automatically if not specified.")
	parser.add_argument('--debug', action='store_true')

	return parser.parse_args()
	
def get_create_time(path):
	t = os.stat(path).st_ctime
	return dt.fromtimestamp(t)

def parse_dates(accno_list):

	if len(accno_list) == 0:
		print("ERROR: At least one accession number is required to query Genbank")
		sys.exit(1)
		
	query = Entrez.efetch(db='nucleotide', id=accno_list, rettype='gb', retmode='text')

	collection_dates = {}

	expr_accno = re.compile('ACCESSION\s+([^\s]+)')
	expr_date = re.compile('=\"(.+)\"')

	# if collection dates are not present, use the oldest submitted journal date as closest estimate
	expr_date_backup = re.compile("Submitted \((.+)\)")

	accno = ''
	backup_date = ''
	for line in query:
		line = line.strip("\n ")

		if line.startswith("ACCESSION"):
			# if the main collection date could not be found and the backup date is available,
			# then use the backup date instead 
			if accno != '' and backup_date != '' and collection_dates[accno] == '':
				collection_dates[accno] = backup_date

			# start of the next gb file 
			print(line)
			accno = expr_accno.search(line).group(1)
			collection_dates[accno] = ''
			backup_date = ''

		# collect and parse backup dates from journal entries
		if line.startswith('JOURNAL   Submitted'):
			date_search = expr_date_backup.search(line).group(1)
			date_search = dt.strptime(date_search, '%d-%b-%Y')

			if isinstance(backup_date, str):
				backup_date = date_search
			else:
				if date_search < backup_date:
					backup_date = date_search

		# retrieve collection dates (preferred date choice)
		if line.startswith("/collection_date"):
			date = expr_date.search(line).group(1)
			collection_dates[accno] = date
	
	if accno != '' and backup_date != '' and collection_dates[accno] == '':
		collection_dates[accno] = backup_date

	# types of genbank date formats
	# written this way with extendability in mind
	# (expression, parse_format, output_format)
	expr_dates = [
		(re.compile('\d{2}-[A-z]+-\d{4}'), '%d-%b-%Y', '%Y-%m-%d'),
		(re.compile('^[A-z]+-\d{4}'), '%b-%Y', '%Y-%m'),
		(re.compile('^\d{4}$'), '%Y', '%Y'),
	]

	for accno, date in collection_dates.items():
		if isinstance(date, dt):
			collection_dates[accno] = date.strftime('%Y-%m-%d')

		elif isinstance(date, str):
			for expr, in_fmt, out_fmt in expr_dates:
				if expr.match(date):
					collection_dates[accno] = dt.strptime(date, in_fmt).strftime(out_fmt)
					break  
		else:
			print(accno)
			print(date)
			print("NOT VALID FORMAT")

	return collection_dates

#%%
def main():
	args = init_parser()

	# parse sequences in FASTA format 
	seqs = list(SeqIO.parse(args.fasta, 'fasta'))

	# retrieve the accession numbers from the FASTA file 
	header_list = pd.Series([x.id for x in seqs])

	# split up the target sequences and the reference sequences
	targets = header_list.str.contains("\|sample")
	background = header_list.str.contains("\|background")

	# create a dictionary of full_header : accession number
	ref_headers = header_list.loc[(~targets)&(~background)]
	ref_headers = {x : x.split(args.header_delim)[args.accno_pos].split(".")[0] for x in ref_headers}

	# extract target sequences
	target_headers = header_list.loc[targets].tolist()
	background_headers = header_list.loc[background].tolist()

	# retrieve the collection dates 
	dates_dict = parse_dates(list(ref_headers.values()))

	# reference dates are collected from Genbank parsing
	ref_dates = [(x, dates_dict[y]) for x, y in ref_headers.items()]

	# target dates are assumed to be current time for simplicity
	target_dates = [(x, dt.today().strftime('%Y-%m-%d')) for x in target_headers]
	background_dates = [(x, x.split("|")[-1]) if x.count("|") == 4 else dt.today().strftime('%Y-%m-%d') for x in background_headers]

	# create two separate dataframes for references + target sequences 
	ref_df = pd.DataFrame(ref_dates, columns=['accno','date'])
	target_df = pd.DataFrame(target_dates, columns=['accno','date'])
	background_df = pd.DataFrame(background_dates, columns=['accno','date'])

	if args.debug:
		print("PARSING ACCNOS:")
		print(list(ref_headers.values()))
		print("PARSING OUTPUT:")
		print(dates_dict)
		print("FINAL OUTPUT:")
		print(ref_df)
		print(target_df)

	if not args.output:
		output_basename = os.path.basename(args.fasta).split(".")[0]
		outpath = output_basename + "_dates.tsv"
	else:
		outpath = args.output

	# output the final combined dataframe
	final_df = pd.concat([ref_df, target_df, background_df]) 
	final_df.to_csv(outpath, sep='\t', index=False, header=False)

if __name__ == '__main__':
	main()