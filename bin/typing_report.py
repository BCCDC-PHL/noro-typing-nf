#!/usr/bin/env python3
#%%
import os, sys
import pandas as pd 
import numpy as np 
import argparse
import itertools

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
	# Read in the BLAST results file and select the mentioned cols
	blast_df = pd.read_csv(path, sep='\t')

	# Checking to parse columns to appropriate datatypes
	float_cols = ['pident', 'bitscore', 'prop_covered', 'composite', 'slen', 'length']
	blast_df[float_cols] = blast_df[float_cols].astype(float).round(1)

	# Sorting the dataframe
	blast_df = blast_df.sort_values(['sample_name','composite'], ascending=False)
	blast_df['sseqid'] = blast_df['sseqid'].str.split("|").str[0]

	# Extracting contig length & depth information from 'qseqid'
	blast_df['contig_len'] = blast_df['qseqid'].str.split("_").str[3].astype(int)
	blast_df['contig_depth'] = blast_df['qseqid'].str.split("_").str[-1].astype(float)

	blast_df['qseqid'] = blast_df['qseqid'].str.split("_length_").str[0]

	# Renaming specific columns for clarity
	blast_df = blast_df.rename({'length': 'hsp_len', 'sseqid':'ref', },axis=1)
	# Final sorting on specified columns before return
	blast_df = blast_df.sort_values(['sample_name','qseqid','composite'])

	blast_df.columns = blast_df.columns.map(lambda x: x if x == 'sample_name' else scheme + "_" + x)

	return blast_df

def merge_blast_dfs(gtype_df, ptype_df):

	# Initialize a set to hold unique 'sample_name'
	sample_list = set(gtype_df['sample_name']).union(ptype_df['sample_name'])

	merged = gtype_df.merge(ptype_df, left_on=['sample_name','g_qseqid'], right_on=['sample_name','p_qseqid'], how='outer')

	merged['combined_score'] = (merged['g_composite'] + merged['p_composite']).round(2)

	# Extract the rows that have both 'g_type' and 'p_type' calls
	full_calls = merged.loc[(~merged['g_type'].isna())&(~merged['p_type'].isna())].copy()

	# Select rows with only one type of call
	partial_gtype = merged.loc[(~merged['g_type'].isna())&(merged['p_type'].isna())].copy()
	partial_ptype = merged.loc[(merged['g_type'].isna())&(~merged['p_type'].isna())].copy()

	# Calculate difference between overall sample list and the ones in 'full calls'
	sample_list = pd.Series(list(sample_list))
	missing_samples = sample_list[~sample_list.isin(full_calls['sample_name'])]

	if len(missing_samples) > 0:

		# Find and keep only the record with max composite score per sample for both type calls
		idxmax = partial_gtype.groupby('sample_name')['g_composite'].idxmax()
		partial_gtype = partial_gtype.loc[idxmax].dropna(axis=1)
		idxmax = partial_ptype.groupby('sample_name')['p_composite'].idxmax()
		partial_ptype = partial_ptype.loc[idxmax].dropna(axis=1)

	# Merge the 'partial' dataframes
	partial_merged = partial_gtype.merge(partial_ptype, on='sample_name', how='outer')
	partial_merged['g_qseqid'] = partial_merged['g_qseqid'].astype(str) + '__' + partial_merged['p_qseqid'].astype(str)

	# Keep only the merged records where sample_name was missing earlier
	partial_merged = partial_merged.loc[partial_merged['sample_name'].isin(missing_samples)]

	full_calls = pd.concat([full_calls, partial_merged])
	
	# Renaming and cleaning up
	full_calls = full_calls.rename({'g_qseqid':'qseqid'},axis=1).drop('p_qseqid',axis=1)

	#full_calls.columns = full_calls.columns.map(lambda x : x.lstrip('p_') if )

	# Combine typing calls into a single 'combined_type' column
	full_calls['g_type'] = full_calls['g_type'].astype(str)
	full_calls['p_type'] = full_calls['p_type'].astype(str)

	#full_calls['combined_type'] = full_calls[['g_type','p_type']].apply(lambda x : "|".join([a+'_'+b for a, b in itertools.zip_longest(x['g_type'].split("|"), x['p_type'].split("|"),fillvalue="NA")]), axis=1)
	full_calls['combined_type'] = full_calls['g_type'] + "_" + full_calls['p_type']
	full_calls = full_calls.drop(['g_type','p_type'], axis=1)

	return full_calls

def sort_df(df, type_col, score_col, id_col='sample_name'):

	df_non_na = df.loc[~df[score_col].isna()]
	df_na = df.loc[df[score_col].isna()]

	# Remove redundancy -- exclude lower quality duplicate typing calls 
	idxmax = df_non_na.groupby([id_col, type_col])[score_col].idxmax()
	df_non_na = df_non_na.loc[idxmax]

	df_na = df_na.loc[~df_na[id_col].isin(df_non_na[id_col])]
	df = pd.concat([df_non_na, df_na])
	# Sort by the appropriate metric
	df = df.sort_values([id_col, score_col],ascending=[True, False])

	return df

def collapse_df(df, id_col):
	# Collapse all rows into one per sample 
	df = df.groupby(id_col).agg(lambda x: '|'.join(map(str, x))).reset_index()
	return df

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
def pass_fail(df, warn_cover=80, fail_cover=50):

	def check(pident, prop_covered):
		if all([np.isnan(pident), np.isnan(prop_covered)]):
			return 'FAIL'
		if any([pident < 80, prop_covered < 80]):
			return 'WARN'
		else:
			return "PASS"
	

	for prefix in ['global_','p_','g_']:
		subset = df[df.columns[df.columns.str.startswith(prefix)]].copy()
		subset = subset.drop(subset.columns[subset.columns.str.endswith(('type','qseqid'))],axis=1)
		subset = subset.apply(max_value)
		newcol = subset.apply(lambda x : check(x[prefix+"pident"], x[prefix+"prop_covered"]), axis=1)
		df.insert(3, prefix+"status", newcol)


#%%
def build_blast_df(gblast_path, pblast_path, global_blast_path):

	gtype_df = parse_blast(gblast_path,'g')
	ptype_df = parse_blast(pblast_path,'p')
	print(gtype_df['g_pident'])

	print(ptype_df['p_pident'])
	global_df = parse_blast(global_blast_path,'global')

	# Remove redundant contig columns 
	gtype_df = gtype_df.drop(gtype_df.columns[gtype_df.columns.str.contains('contig')], axis=1)
	ptype_df = ptype_df.drop(ptype_df.columns[ptype_df.columns.str.contains('contig')], axis=1)

	# merge the individual typing calls into one 
	combined_df = merge_blast_dfs(gtype_df, ptype_df)

	# sort dataframes 
	global_df = sort_df(global_df, 'global_type', 'global_composite')
	combined_df = sort_df(combined_df, 'combined_type', 'combined_score')

	# collapse dataframes to show the top 3 results
	global_df = collapse_df(global_df, 'sample_name')
	combined_df = collapse_df(combined_df, 'sample_name')

	full_blast = global_df.merge(combined_df, on='sample_name', how='left')
	full_blast.columns = full_blast.columns.map(lambda x : x.lstrip('global_') if 'contig' in x else x)
	full_blast = full_blast.rename({'qseqid':'combined_qseqid'})


	full_blast = full_blast.drop(full_blast.columns[full_blast.columns.str.endswith(('ref', 'qseqid'))], axis=1)

	full_blast.insert(1, 'combined_type', full_blast.pop('combined_type'))
	full_blast.insert(2, 'global_type', full_blast.pop('global_type'))
	
	str_cols = full_blast.columns[full_blast.dtypes == 'object']
	float_cols = full_blast.columns[full_blast.dtypes == 'float']

	full_blast[str_cols].fillna('NA', inplace=True)
	full_blast[float_cols].fillna(np.nan, inplace=True)

	pass_fail(full_blast)

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