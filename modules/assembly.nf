process assembly {

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

process run_shovill {
	tag {sample_id}

    publishDir "${params.outdir}/assembly", pattern: "${sample_id}.contigs.fa", mode:'copy'
	publishDir "${params.outdir}/assembly/full", pattern: "${sample_id}.spades.tar.gz", mode:'copy'

	input:
	tuple val(sample_id), path(reads_1), path(reads_2)

	output: 
	tuple val(sample_id), path("${sample_id}.contigs.fa")
	
	"""
	shovill --outdir --cpus ${task.cpus} --R1 ${reads_1} --R2 ${reads_2}
	"""

}

process run_metabat {

	errorStrategy 'ignore'

	tag {sample_id}

	conda "metabat2"

	publishDir "${params.outdir}/assembly/binning/mb", pattern: "${sample_id}*", mode:'copy'

	input:
	tuple val(sample_id), path(contig_fasta), path(aligned_bam), path(bam_index)

	output: 
	tuple val(sample_id), path("${sample_id}*")
	
	"""
	runMetaBat.sh ${contig_fasta} ${aligned_bam}
	mkdir -p metabat
	mv ${sample_id} ${sample_id}_metabat
	"""
}

process run_concoct {
	tag {sample_id}

	errorStrategy 'ignore'

	label 'medium'

	conda "${projectDir}/environments/concoct2.yaml"

	publishDir "${params.outdir}/assembly/binning", pattern: "${sample_id}_bins", mode:'copy'
	publishDir "${params.outdir}/assembly/binning/raw", pattern: "${sample_id}.concoct.tar.gz", mode:'copy'

	input:
	tuple val(sample_id), path(contig_fasta), path(aligned_bam), path(bam_index)

	output: 
	tuple val(sample_id), path("${sample_id}_bins"), path("${sample_id}.concoct.tar.gz")
	
	"""
	mkdir -p ${sample_id}
	mkdir -p ${sample_id}_bins

	cut_up_fasta.py ${contig_fasta} -c 10000 -o 0 --merge_last -b ${sample_id}/${sample_id}_contigs_10K.bed > ${sample_id}/${sample_id}_contigs_10K.fa
	
	# Generate coverage depth 
	concoct_coverage_table.py ${sample_id}/${sample_id}_contigs_10K.bed ${aligned_bam} > ${sample_id}/${sample_id}_coverage_table.tsv
	
	# Execute CONCOCT
	concoct --threads ${task.cpus} --composition_file ${sample_id}/${sample_id}_contigs_10K.fa --coverage_file ${sample_id}/${sample_id}_coverage_table.tsv -b ${sample_id}/${sample_id}
	
	# Merge sub-contig clustering into original contig clustering
	merge_cutup_clustering.py ${sample_id}/${sample_id}_clustering_gt1000.csv > ${sample_id}/${sample_id}_clustering_merged.csv
	
	# Parse bins into different files
	extract_fasta_bins.py ${contig_fasta} ${sample_id}/${sample_id}_clustering_merged.csv --output_path ${sample_id}_bins

	tar -czvf ${sample_id}.concoct.tar.gz ${sample_id}
	"""
}