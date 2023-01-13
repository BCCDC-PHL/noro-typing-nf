# %%
import os
import pandas as pd
from Bio import SeqIO
from Bio import Entrez

# %%
Entrez.email = ''

# %%
base = ''
os.chdir(base)
print(os.getcwd())

# %%
df = pd.read_csv('ptypes.csv')
df['combined'] = df['accession'] + '|' + df['p-type'] + '|' + df['strain']

# %% fetch entries from NCBI Genbank
result = Entrez.efetch(db='nucleotide', id=df['accession'].to_list(), rettype='fasta',retmode='text')



#%% write output to FASTA
with open('entrez.fasta', 'w') as outfile:
	seqs = []
	for i in result:

		outfile.write(i)


seqs = SeqIO.to_dict(SeqIO.parse('entrez.fasta','fasta'))

# %%
df['accession'].duplicated().sum()

#%% 
df.loc[df['accession'].duplicated()]
df = df.drop_duplicates('accession')

#%% setup dictionary with key = old header, value = new header
accno_dict = df.set_index('accession')['combined'].to_dict()

#%% clean headers using accno_dict 
for header, seq in seqs.items():
	header = header.split(".")[0]
	seq.name = '' 
	seq.id = accno_dict[header]
	print(header)
	seq.description = '' 

# %% write final database to FASTA 
SeqIO.write(list(seqs.values()), 'norovirus-ptype-db.fasta', 'fasta')
