process assembly {

	errorStrategy 'ignore'

	tag {sample_id}

    publishDir "${params.outdir}/assembly", pattern: "${sample_id}.contigs.fa", mode:'copy'
	publishDir "${params.outdir}/assembly/full", pattern: "${sample_id}.spades.tar.gz", mode:'copy'

	input:
	tuple val(sample_id), path(reads_1), path(reads_2)

	output: 
	tuple val(sample_id), path("${sample_id}.contigs.fa")
	
	"""
	rnaviralspades.py -1 ${reads_1} -2 ${reads_2} -o ${sample_id} &&
	cp ${sample_id}/contigs.fasta ./${sample_id}.contigs.fa &&
	tar -czvf ${sample_id}.spades.tar.gz ${sample_id}/*
	"""

}

process run_metabat {

	tag {sample_id}

	publishDir "${params.outdir}/assembly/binning", pattern: "${sample_id}*", mode:'copy'

	input:
	tuple val(sample_id), path(contig_fasta), path(aligned_bam)

	output: 
	tuple val(sample_id), path("${sample_id}.metabat.tar.gz")
	
	"""
	runMetaBat.sh ${contig_fasta} ${aligned_bam}
	mkdir -p metabat
	mv ${sample_id} metabat
	tar -czvf ${sample_id}.metabat.tar.gz metabat/*
	"""
}

process run_concoct {
	tag {sample_id}

	publishDir "${params.outdir}/assembly/binning", pattern: "${sample_id}.spades.tar.gz", mode:'copy'

	input:
	tuple val(sample_id), path(contig_fasta), path(aligned_bam)

	output: 
	tuple val(sample_id), path("${sample_id}.contigs.fa")
	
	"""
	cut_up_fasta.py ${contig_fasta} -c 10000 -o 0 --merge_last -b ${sample_id}_contigs_10K.bed > cut_fasta/${sample_id}_contigs_10K.fa
	
	# Generate coverage depth 
	concoct_coverage_table.py cut_fasta/${sample_id}_contigs_10K.bed ${aligned_bam} > coverage/${sample_id}_coverage_table.tsv
	
	# Execute CONCOCT
	concoct --threads ${task.cpus} --composition_file cut_fasta/${sample_id}_contigs_10K.fa --coverage_file coverage/${sample_id}_coverage_table.tsv -b output/${sample_id}
	
	# Merge sub-contig clustering into original contig clustering
	merge_cutup_clustering.py output/${sample_id}_clustering_gt1000.csv > output/${sample_id}_clustering_merged.csv
	
	# Parse bins into different files
	extract_fasta_bins.py ${contig_fasta} ${sample_id}_clustering_merged.csv --output_path ${sample_id}
	"""
}