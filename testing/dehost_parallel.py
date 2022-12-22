
#%%
import multiprocessing as mp
import pysam 
import os 
from collections import defaultdict

def read_alignment(path):
	if not path:
		infile = pysam.AlignmentFile('-', 'r')
	elif path.endswith('.sam'):
		infile = pysam.AlignmentFile(path, 'r')
	else:
		infile = pysam.AlignmentFile(path, 'rb')
	return infile
	
def write_alignment(bam_output_path, header, mode='w'):
	"""
	Function that returns a handle that writes to a new output sam file
	"""
	if bam_output_path:
		output_bam = pysam.AlignmentFile(bam_output_path, mode, header=header)
	else:
		output_bam = pysam.AlignmentFile('-', mode, header=header)
	return output_bam


#%%
def split_sam(filepath, threads):

	input_sam = read_alignment(filepath)
	header = input_sam.header
	N = input_sam.count()

	input_sam.close()

	part_size = N // threads + 1

	collect = []
	read_generators = []
	read_counts = []
	count = 0 

	output_names = [os.path.join(os.path.dirname(filepath), '.'+os.path.basename(filepath)+'.'+str(i)) for i in range(threads)]

	current_outfile = None

	input_sam = read_alignment(filepath)
	
	# sequentially dump reads into multiple generators 
	for n, read in enumerate(input_sam.fetch()):
		if n % part_size == 0:
			if current_outfile:
				current_outfile.close()
				read_counts.append(count)
				count = 0
		
			current_outfile = write_alignment(output_names[n // part_size], input_sam.header)

		count += 1
		current_outfile.write(read)
	
	read_counts.append(count)
	current_outfile.close()

	# count += 1
	# collect.append(read.to_string())

	# if n % part_size == 0: 
	# 	#read_generators.append((x for x in collect))
	# 	read_generators.append( collect)
	# 	read_counts.append(count)

	# 	# reset
	# 	count = 0
	# 	collect = []

	# append the final generator and count 
	# read_generators.append((x for x in collect))
	# read_generators.append( collect)
	# read_counts.append(count)
	
	return output_names, read_counts
#%%
def filter_reads(input_path, output_path, viral_ref_names):
	# if read isn't mapped or mapped to viral reference contig name

	human_reads = 0
	unmapped_reads = 0
	viral_dict = dict(zip(viral_ref_names, [0] * len(viral_ref_names)))

	input_bam = read_alignment(input_path)

	output_bam = write_alignment(output_path, input_bam.header)

	print('starting async filtering')
	# iterate over input from BWA
	for n, read in enumerate(input_bam.fetch()):
		# only look at primary alignments
		if not read.is_supplementary and not read.is_secondary:
			if read.reference_name in viral_dict:
				output_bam.write(read)
				viral_dict[read.reference_name] += 1
			elif read.is_unmapped:
				output_bam.write(read)
				unmapped_reads +=1
			else:
				human_reads += 1

	output_bam.close()

	print("finished filtering")
	viral_reads = sum(viral_dict.values())
	total_reads = viral_reads + human_reads + unmapped_reads

	if total_reads > 0:
		print(f"viral read count = {viral_reads} "
			f"({viral_reads/total_reads * 100:.2f}%)", file=sys.stderr)
		print(viral_dict, file=sys.stderr)
		print(f"human read count = {human_reads} "
			f"({human_reads/total_reads * 100:.2f}%) ", file=sys.stderr)
		print(f"unmapped read count = {unmapped_reads} "
			f"({unmapped_reads/total_reads * 100:.2f}%)", file=sys.stderr)
	else:
		print(f"ERROR: No reads remaining", file=sys.stderr)
	
	return viral_dict, human_reads, unmapped_reads

# Function to filter human reads
# Derived from ncov2019-artic-nf 
#%%
def filter_reads_parallel(input_sam_fp, output_bam_fp, viral_ref_names, threads):
	
	print('beginning to split sam')
	
	tmp_files, counts = split_sam(input_sam_fp, threads=threads)

	print('running parallel')

	# results = []
	# def log_results(res):
	# 	global results
	# 	results.append(res)

	# for n, i in enumerate(generator_list[0]):
	# 	print(i.reference_name)
	# 	if n == 100: 
	# 		break

	output_files = [x.replace('.sam','.filter.sam') for x in tmp_files]

	jobs = []
	with mp.Pool(threads) as pool:

		for n, (infile, outfile) in enumerate(zip(tmp_files, output_files)):
			print('process ', n)
			job = pool.apply_async(filter_reads, args=(infile, outfile, viral_ref_names))
			jobs.append(job)

		for job in jobs:
			print(type(job))
			job.get()

	for tmp in tmp_files:
		os.remove(tmp)

	print("waiting...")	

	return 'success'
