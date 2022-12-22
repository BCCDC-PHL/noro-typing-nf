# %%
import os
import pandas as pd
from Bio import SeqIO
from Bio import Entrez

# %%
base = '/home/john.palmer/work/norovirus/blast'
os.chdir(base)

# %%
df = pd.read_csv('subtypes.csv')
df['combined'] = df['accession'] + '|' + df['genotype'] + '|' + df['strain']

# %%
Entrez.email = 'john.palmer1288@gmail.com'
Entrez.api_key = os.environ['NCBI_API_KEY']
# %%s
os.getcwd()

# %%
result = Entrez.efetch(db='nucleotide', id=df['accession'].to_list(), rettype='fasta',retmode='text')
with open('entrez.fasta', 'w') as outfile:
	for i in result:
		outfile.write(i)
seqs = SeqIO.to_dict(SeqIO.parse('entrez.fasta','fasta'))

# %%
df['accession'].duplicated().sum()
df.loc[df['accession'].duplicated()]
df = df.drop_duplicates('accession')

# %%
list(seqs.values())[0:5]

# %%
# %%
for header, seq in seqs.items():
	header = header.split(".")[0]
	seq.name = '' 
	seq.id = df.loc[df['accession'] == header, 'combined'].squeeze()
	print(header)
	seq.description = '' 

# %%
SeqIO.write(list(seqs.values()), 'norovirus-db-v2.fasta', 'fasta')


