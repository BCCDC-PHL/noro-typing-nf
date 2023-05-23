import os, sys
import pandas as pd 
from glob import glob 

dflist = []

for file in glob(os.path.join(sys.argv[1], '*report')):
	sample_name = os.path.basename(file).split(".")[0]
	print(sample_name)
	df = pd.read_csv(file, sep='\t', header=None)
	df.columns = ['proportion', 'subcount', 'count', 'abbrev', 'id', 'species']

	df['species'] = df['species'].str.strip()

	df = df.loc[df['count'] > 200]
	df.insert(0, 'sample', sample_name)

	target = df['species'].str.startswith(('Norwalk','Norovirus'))
	total = df.loc[target, 'count'].astype(int).sum()
	df = df.drop(df.index[target], axis=0).reset_index(drop=True)

	main = df['species'].str.startswith(('Norovirus', 'unclassified','Homo sapiens'))
	df.loc[~main, 'species'] = 'other'

	df.loc[df.shape[0]] = [sample_name, -1, -1, total, 'NA', -1, 'Norovirus']
	

	dflist.append(df)

final = pd.concat(dflist)
final = final.reset_index(drop=True)
print(final.shape)
print(final.head())

final.to_csv(sys.argv[2], index=False)