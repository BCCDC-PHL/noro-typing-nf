#!/usr/bin/env python3
#%%
import os, sys
import pandas as pd 
import numpy as np 
import argparse
import itertools
import re

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_colwidth', 500)
#%%
def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--sample-list', required=True, help='Full list of sample names')
	parser.add_argument('-a', '--quast', required=False, help='Path to the Quast assembly QC report TSV')
	parser.add_argument('-q', '--qc', required=True, help='Path to the custom QC files')
	parser.add_argument('-g', '--gblast',  required=True, help='Path to genotype BLAST results')
	parser.add_argument('-p', '--pblast', required=True, help='Path to ptype BLAST results')
	parser.add_argument('-f', '--globalblast', required=False, help='Path to full BLAST results')
	parser.add_argument('-o', '--outpath', required=True, help='Name and path of final output report')

	return parser

def parse_sample_file(path):
	with open(path, 'r') as infile:
		samples = [x.strip() for x in infile.readlines()]

	samples_df = pd.DataFrame({'sample_name':samples})

	return samples_df

def parse_qc(path):
	qc_df = pd.read_csv(path)
	return qc_df

def parse_quast(path):
	quast_df = pd.read_csv(path, sep='\t')
	quast_df = quast_df.transpose()
	quast_df.columns = quast_df.iloc[0,:]
	quast_df = quast_df.drop(quast_df.index[0])
	drop_list = ['# contigs (>= 25000 bp)', '# contigs (>= 10000 bp)', '# contigs (>= 50000 bp)', "# N's per 100 kbp",'GC (%)'
		  		#'Total length (>= 10000 bp)','Total length (>= 25000 bp)','Total length (>= 50000 bp)', 
			]
	quast_df = quast_df.drop(drop_list, axis=1)
	quast_df = quast_df.drop(quast_df.columns[quast_df.columns.str.contains('Total length \(')],axis=1)
	quast_df = quast_df.reset_index(names='sample_name')
	quast_df['sample_name'] = quast_df['sample_name'].str.split(".").str[0]
	# quast_df['Assembly'] = quast_df['Assembly'].str.split(".").str[0]

	return quast_df

#%%
def parse_blast(path, scheme):
	# Define column headers
	cols = 'sample_name qseqid sseqid slen pident length prop_covered type bitscore composite'.split()
	# Read in the BLAST results file and select the mentioned cols
	blast_df = pd.read_csv(path, sep='\t')[cols]

	# Sorting the dataframe
	blast_df = blast_df.sort_values(['sample_name','composite'], ascending=False)

	# Extracting contig length & depth information from 'qseqid'
	blast_df['contig_len'] = blast_df['qseqid'].str.split("_").str[3]
	blast_df['contig_depth'] = blast_df['qseqid'].str.split("_").str[-1]

	# Checking to parse columns to appropriate datatypes
	int_cols = ['contig_len', 'slen', 'length']
	blast_df[int_cols] = blast_df[int_cols].astype(int)
	float_cols = ['contig_depth', 'prop_covered', 'composite']
	blast_df[float_cols] = blast_df[float_cols].astype(float).round(2)

	# Clean sequence IDs
	blast_df['sseqid'] = blast_df['sseqid'].str.split("|").str[0]
	blast_df['qseqid'] = blast_df['qseqid'].str.split("_length_").str[0]

	# Renaming specific columns for clarity
	blast_df = blast_df.rename({'length': 'hsp_length','sseqid':'reference', },axis=1)

	# Final sorting on specified columns before return
	blast_df = blast_df.sort_values(['sample_name','qseqid','composite']).reset_index(drop=True)

	blast_df.columns = blast_df.columns.map(lambda x: x if x == 'sample_name' else scheme + "_" + x)

	return blast_df


def merge_blast(dfs):

	# Initialize a set to hold unique 'sample_name'
	sample_list = set()

	# Loop over the two dataframes ('g' and 'p') in the list
	for n, scheme in enumerate(['g','p']):
		blast_df = dfs[n]
		sample_list = sample_list.union(set(blast_df['sample_name']))

		# For the first iteration, copy the dataframe
		if n == 0:
			merged = blast_df.copy()
		else:
			# For the second iteration, merge the dataframes
			merged = merged.merge(blast_df, left_on=['sample_name','g_qseqid'], right_on=['sample_name','p_qseqid'], how='outer')

	# Renaming and cleaning up
	merged = merged.sort_values(['sample_name','g_qseqid'])
	
	merged = merged.reset_index(drop=True)

	return merged

def handle_partials(df, score_col):
	# Extract the rows that have both 'g_type' and 'p_type' calls
	complete = df.loc[(~df['g_type'].isna()) & (~df['p_type'].isna())].copy()
	complete.insert(1, 'status', ['COMPLETE'] * complete.shape[0])
	

	# Select rows with only one type of call
	partial_gtype = df.loc[(df['p_type'].isna())].copy()
	partial_ptype = df.loc[(df['g_type'].isna())].copy()

	# Find and keep only the record with max composite score per sample for both type calls
	idxmax1 = partial_gtype.groupby('sample_name')['g_'+score_col].idxmax()
	partial_gtype = partial_gtype.loc[idxmax1].dropna(axis=1)
	idxmax2 = partial_ptype.groupby('sample_name')['p_'+score_col].idxmax()
	partial_ptype = partial_ptype.loc[idxmax2].dropna(axis=1)

	# Merge the 'partial' dataframes
	partial_merged = partial_gtype.merge(partial_ptype, on='sample_name', how='inner')
	#partial_merged['g_qseqid'] = partial_merged['g_qseqid'].astype(str) + '__' + partial_merged['p_qseqid'].astype(str)
	partial_merged.insert(1, 'status', ['COMBINED'] * partial_merged.shape[0])

	combined_samples = set(partial_gtype.sample_name).intersection(partial_ptype.sample_name)
	merged_rows = set(partial_gtype.index[partial_gtype.sample_name.isin(combined_samples)]).union(set(partial_ptype.index[partial_ptype.sample_name.isin(combined_samples)]))

	remainder_rows = list(set(df.index) - set(complete.index) - merged_rows)
	remainder = df.loc[remainder_rows]
	remainder.insert(1, 'status', ['PARTIAL'] * remainder.shape[0])

	cleaned = pd.concat([complete, partial_merged, remainder]).sort_values(['sample_name','g_qseqid']).reset_index(drop=True)
	cleaned.insert(2, "synth_score", (cleaned['g_'+score_col].fillna(0) + cleaned['p_'+score_col].fillna(0)).round(2))

	return cleaned

def zip_columns(df, columns, delim='|', join_char='__'):
	df[columns] = df[columns].fillna("NA")
	return df[columns].apply(lambda x : delim.join([a+join_char+b for a, b in itertools.zip_longest(x[columns[0]].split(delim), x[columns[1]].split(delim),fillvalue="NA")]), axis=1)


def clean_dataframe(df,  score_col='bitscore'):

	df = handle_partials(df, score_col)

	df[['g_type', 'p_type']] = df[['g_type', 'p_type']] .fillna("NA")

	idxmax = df.groupby(['sample_name','g_type', 'p_type'])['synth_score'].idxmax()
	df = df.loc[idxmax]

	# Sort by the appropriate metric
	df = df.sort_values(['sample_name','synth_score'],ascending=[True, False])
	df.insert(3, "contigs", zip_columns(df, ['g_qseqid','p_qseqid']))
	
	condensed = df.groupby('sample_name').agg(lambda x: '|'.join(map(str, x))).reset_index().reset_index(drop=True).copy()
	df.insert(3, "synth_type", df['g_type'] + '__' + df['p_type'])

	condensed.insert(3, "synth_type", zip_columns(df, ['g_type','p_type']))

				  
	
	#df = df.drop(['g_type','p_type'], axis=1)
	#condensed = condensed.drop(['g_type','p_type'], axis=1)

	single = condensed.apply(lambda x : x if x.dtype != 'object' else x.str.split("\|").str[0])

	return df, condensed


#%%
def max_value(string, delim='|'):
	if not isinstance(string, pd.Series):
		print("ERROR: Incorrect input when parsing top hits.")
		exit(1)

	if string.dtype != 'object':
		print("ERROR: Incorrect series data type")
		exit(0)

	string = string.str.split(delim).str[0]

	try: 
		string = string.astype(float)
	except ValueError as e:
		pass
	return string

#%%
def build_blast_df(gblast_path, pblast_path, global_blast_path):

	gtype_df = parse_blast(gblast_path,'g')
	ptype_df = parse_blast(pblast_path,'p')

	# merge the individual typing calls into one 
	merged = merge_blast([gtype_df, ptype_df])	
	complete, condensed = clean_dataframe(merged, score_col='bitscore')



	if global_blast_path:
		global_df = parse_blast(global_blast_path,'global')
		# sort and collapse dataframes to show the top 3 results
		global_df = sort_collapse_df(global_df, 'global_type', 'global_composite')

		full_blast = global_df.merge(synth_df, on='sample_name', how='left')
		full_blast.insert(1, 'global_type', full_blast.pop('global_type'))

	else:
		full_blast = synth_df


	# contig_cols = {}
	# for col in full_blast.columns:
	# 	if "contig" in col:
	# 		search = re.search("^[^_]+_(.+)", col)
	# 		simple_name = search.group(1) if search else None
	# 		if simple_name and simple_name not in contig_cols:
	# 			contig_cols[simple_name] = col
	# contig_cols = {y:x for x, y in contig_cols.items()}

	# full_blast = full_blast.rename(contig_cols, axis=1)
	# full_blast = full_blast.drop(full_blast.columns[full_blast.columns.str.contains("_contig_")], axis=1)


	full_blast = full_blast.fillna('NA')

	full_blast = full_blast.drop(full_blast.columns[full_blast.columns.str.endswith(('ref', 'qseqid'))], axis=1)
	full_blast.insert(1, 'synth_type', full_blast.pop('synth_type'))

	top_blast = full_blast.copy()
	top_blast = top_blast.apply(lambda x : x.str.split("|").str[0] if x.dtype == 'object' else x)

	return full_blast, top_blast

#%%
def main(args):
	# Parse the custom QC report containing mapping depth and coverage information 
	qc_df = parse_qc(args.qc)

	# Build the BLAST typing dataframe 
	top_three_df, top_df = build_blast_df(args.gblast, args.pblast, args.globalblast)

	# Load in the samples from the sample list txt file 
	sample_df = parse_sample_file(args.sample_list)

	# Build the top 3 report
	top_three_df = sample_df.merge(top_three_df, on='sample_name', how='left')
	top_three_df = top_three_df.merge(qc_df, on='sample_name', how='left')
	top_three_df.to_csv(args.outpath.split(".")[0] + "_top3.tsv", sep='\t', index=False)

	# Build the top 1 report 
	top_df = sample_df.merge(top_df, on='sample_name', how='left')
	top_df = top_df.merge(qc_df, on='sample_name', how='left')
	top_df.to_csv(args.outpath, sep='\t', index=False)

#%%
if __name__ == '__main__':
	parser = init_parser()
	main(parser.parse_args())