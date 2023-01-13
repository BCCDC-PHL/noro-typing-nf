import os, sys, re
from Bio import SeqIO
from Bio import Entrez
from collections import defaultdict, Counter

Entrez.email = 'john.palmer1288@gmail.com'
Entrez.api_key = os.environ['NCBI_API_KEY']


def reformat_headers(seq_list):
	expr = re.compile('\|(.+\.)\d_prot.+_(\d)')

	final_seqs = seq_list.copy()

	for seq in final_seqs:
		search = expr.search(seq.id)

		num_groups = sum(x is not None for x in search.groups())

		if num_groups != 2:
			print(f"ERROR: {seq.id} failed to be parsed by regex")
			continue 
		
		new_name = "".join(search.groups())

		seq.id = new_name if new_name else seq.id
		seq.name = new_name.split(".")[0] if new_name else seq.name
		seq.description = new_name if new_name else seq.description

	final_seqs = sorted(final_seqs, key=lambda x: x.id)

	return final_seqs


def remove_partial(seq_list, num_cds):
	"""
	Acts on the basis that seq.name is the common accession number of the sequence WITHOUT any decimal ending (MG495082)
	and are common between CDSs from the same GenBank entry (as they are used in filterin)
	"""
	count = Counter((x.name for x in seq_list))

	full_ids = {x for x,y  in count.items() if y == num_cds}

	print("Original Sequence Count: ", len(seq_list))
	print("Partial Sequence Count: ", len(seq_list) - len(full_ids)*num_cds)
	print("Full Sequences Only: ", len(full_ids)*num_cds)

	return [x for x in seq_list if x.name in full_ids]



def collect_cds_regions(seq_list, num_cds):

	seq_dict = defaultdict(dict)

	for record in seq_list:
		seq_dict[record.name][record.id.split(".")[-1]] = record

	final_dict = {}

	for header in seq_dict.keys():
		for idx in range(1, num_cds+1):
			if header not in final_dict:
				final_dict[header] = seq_dict[header][str(idx)]
			else:
				final_dict[header].seq += seq_dict[header][str(idx)].seq
	
	return final_dict

if __name__ == '__main__':
	seqs = list(SeqIO.parse('efetch.fasta','fasta'))

	renamed_seqs = reformat_headers(seqs)

	full_seqs = remove_partial(renamed_seqs, 3)

	final_seqs = collect_cds_regions(full_seqs, 3)

	SeqIO.write(list(final_seqs.values()), 'protein_seqs.fasta','fasta')


