#%%
import os, sys
import pandas as pd 
import numpy as np 
from glob import glob 
import argparse
import itertools

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_colwidth', 500)

def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--names', required=True, help='Full list of sample names')
	parser.add_argument('-a', '--quast', required=True, help='Path to the Quast assembly QC report TSV')
	parser.add_argument('-q', '--qc', required=True, help='Path to all QC TSV files')
	parser.add_argument('-g', '--gblast', nargs='+', required=True, help='All filtered genotype BLAST results')
	parser.add_argument('-p', '--pblast', required=True, help='Path to all ptype BLAST results')

	return parser

def parse_sample_file(path):
	with open(path, 'r') as infile:
		samples = [x.strip() for x in infile.readlines()]
	return samples

def parse_qc(path):
	linelist = []
	first = True

	fullpath = os.path.abspath(path)

	for path in glob(fullpath + "/*qc.csv"):
		with open(path, 'r') as infile:
			lines = [x.strip().split(",") for x in infile.readlines()]
			if first: 
				linelist.append(lines[0])
				first = False
			linelist.append(lines[1])
	
	qc_df = pd.DataFrame(linelist[1:], columns=linelist[0])

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
def parse_blast(paths, mode):

	linelist = []
	first = True

	for path in paths:

		fullpath = os.path.abspath(path)

		with open(fullpath, 'r') as infile:
			lines = [x.strip().split("\t") for x in infile.readlines()]
			if first: 
				linelist.append(lines[0])
				first = False
			linelist += lines[1:]

	cols = 'sample_name qseqid sseqid slen pident length prop_covered cov type composite'.split()
	blast_df = pd.DataFrame(linelist[1:], columns=linelist[0])[cols]

	int_cols = ['slen', 'length']
	blast_df[int_cols] = blast_df[int_cols].astype(int)
	float_cols = ['cov', 'prop_covered', 'composite']
	blast_df[float_cols] = blast_df[float_cols].astype(float).round(2)

	blast_df = blast_df.sort_values(['sample_name','composite'], ascending=False)
	blast_df['sseqid'] = blast_df['sseqid'].str.split("|").str[0]

	max_hits = blast_df['sample_name'].value_counts().max()
	max_hits = 3 if max_hits > 3 else max_hits

	blast_df = blast_df.rename({'length': 'hsp_len', 'cov': 'contig_depth','sseqid':'ref'},axis=1)
	blast_df['contig_len'] = blast_df['qseqid'].str.split("_").str[3].astype(int)

	blast_df = blast_df.groupby('sample_name').agg(lambda x: '|'.join(map(str, x))).reset_index()

	blast_df = blast_df.drop(['hsp_len', 'slen', 'qseqid'], axis=1)
	#print(blast_df.columns)
	blast_df = blast_df[['sample_name', 'ref', 'pident', 'prop_covered', 'type', 'composite', 'contig_depth', 'contig_len']]
	blast_df.columns = blast_df.columns.map(lambda x: mode+"_"+x if x != 'sample_name' else x)

	return blast_df

import numpy as np
#%%
def max_value(string, delim='|'):
	string = string.str.split(delim).str[0]
	string = string.replace({'NA':np.nan})
	return string.astype(float)

def pass_fail(df, warn_cover=80, fail_cover=30):
	cover_cols = ['full_prop_covered', 'g_prop_covered', 'p_prop_covered']   # using lists to make order unambiguous
	hsp_len = ['full_hsp_len', 'g_hsp_len', 'p_hsp_len']  


	bins = [0, fail_cover, warn_cover, 100]

	subset = df[cover_cols].copy()
	subset = subset.apply(max_value)
	subset = subset.apply(lambda x : pd.cut(x, bins, labels=['FAIL','WARN','PASS']))
	
	def classify(full, gtype, ptype):
		if full == 'PASS' or (gtype =='PASS' and ptype =='PASS'):
			return 'PASS'
		elif full != 'FAIL' or (gtype !='FAIL' and ptype != 'FAIL'):
			return 'WARN'
		else:
			return 'FAIL'
	
	final = subset.apply(lambda x : classify(*[x[i] for i in cover_cols]), axis=1)

	return final

def pass_fail2(df, warn_cover=80, fail_cover=30):

	def check(pident, prop_covered, contig_len):
		if all([np.isnan(pident), np.isnan(prop_covered), np.isnan(contig_len)]):
			return 'FAIL'
		if any([pident < 80, prop_covered < 80, contig_len < 1000]):
			return 'WARN'
		else:
			return "PASS"
	

	for prefix in ['full_','p_','g_']:
		subset = df[df.columns[df.columns.str.startswith(prefix)]].copy()
		subset = subset.drop(subset.columns[subset.columns.str.endswith('type')],axis=1)
		subset = subset.apply(max_value)

		newcol = subset.apply(lambda x : check(x[prefix+"pident"], x[prefix+"prop_covered"],x[prefix+"contig_len"]), axis=1)
		df.insert(3, prefix+"status", newcol)


#%%
def merge_blast():
	gblast = glob("/home/john.palmer/work/norovirus/results-may1-4/blastn/gtype/filtered/*tsv")
	pblast = glob("/home/john.palmer/work/norovirus/results-may1-4/blastn/ptype/filtered/*tsv")
	fullblast = glob("/home/john.palmer/work/norovirus/results-may1-4/blastn/global/filtered/*tsv")
	quast_path = '/home/john.palmer/work/norovirus/results-may1-4/qc/assembly/reports/quast_results/report.tsv'
	quast_path = '/home/john.palmer/work/norovirus/results-may1-4/qc/custom/reports/quast_results/report.tsv'

	gtypedf = parse_blast(gblast, 'g')
	ptypedf = parse_blast(pblast, 'p')
	fulldf = parse_blast(fullblast, 'full')


	fullblast = fulldf.merge(gtypedf, on='sample_name', how='left')
	fullblast = fullblast.merge(ptypedf, on='sample_name', how='left')
	fullblast = fullblast.fillna('NA')
	#print(fullblast)
	fullblast['synth_type'] = fullblast[['g_type','p_type']].apply(lambda x : "|".join([a+'_'+b for a, b in itertools.zip_longest(x['g_type'].split("|"), x['p_type'].split("|"),fillvalue="NA")]), axis=1)
	fullblast = fullblast.drop(['g_type','p_type'], axis=1)
	fullblast = fullblast.drop(fullblast.columns[fullblast.columns.str.endswith('ref')], axis=1)

	fullblast.insert(1, 'synth_type', fullblast.pop('synth_type'))
	fullblast.insert(2, 'full_type', fullblast.pop('full_type'))

	pass_fail2(fullblast)
	# fullblast.insert(3, 'status', pass_fail(fullblast))

	return fullblast# , parse_quast(quast_path), parse_qc(qc_path)

#%%
def main(args):
	samples = parse_sample_file(args.names)

	qc_df = parse_qc(args.qc)

	quast_df = parse_quast(args.quast)

	full = merge_blast()

	samples = parse_sample_file('/home/john.palmer/work/norovirus/results-may1-3/samples.txt')

	samples_df = pd.DataFrame({'sample_name':samples})

	fulldf = samples_df.merge(full, on='sample_name', how='left')

	print(fulldf)

	# gtypedf = parse_blast(args.gblast, 'g')
	# ptypedf = parse_blast(args.gblast, 'p')
	# fulldf = parse_blast(args.global, 'full')

#%%
if __name__ == '__main__':
	parser = init_parser()
	main(parser.parse_args())