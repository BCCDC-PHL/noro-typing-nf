process get_coverage { 
	errorStrategy 'ignore'

	tag { sample_id }

	publishDir "${params.outdir}/qc/coverage/bed", mode:'copy'

	input:
	tuple val(sample_id), path(bamfile), path(bam_index)

	output:
	tuple val(sample_id), path("*bed"),emit: coverage_file
	tuple val(sample_id), path("${sample_id}_samtools_provenance.yml"), emit: provenance

	"""
	printf -- "- process_name: samtools\\n" > ${sample_id}_samtools_provenance.yml
	printf -- "  tool_name: samtools\\n tool_version: \$(samtools | sed -n '1 p'))\\n" >> ${sample_id}_samtools_provenance.yml

	samtools depth -a ${bamfile} > ${sample_id}_coverage.bed
	"""
}

process plot_coverage { 
	errorStrategy 'ignore'

	tag { sample_id }

	conda "${projectDir}/environments/plot.yaml"

	publishDir "${params.outdir}/qc/coverage/plots", mode:'copy'

	input:
	tuple val(sample_id), path(coverage_bed)

	output:
	tuple val(sample_id), path("${sample_id}*.pdf")

	"""
	plotting.R ${coverage_bed} ${sample_id}_coverage.pdf
	"""
}

process make_pileup {
	errorStrategy 'ignore'

	tag { sample_id }

	publishDir "${params.outdir}/qc/pileups/ambiguous", pattern: "*ambiguous.tsv", mode:'copy'
	publishDir "${params.outdir}/qc/pileups/full",  pattern: "*pileup.tsv", mode:'copy'

	input:
	tuple val(sample_id), path(reference), path(bam), path(bam_index)

	output:
	tuple val(sample_id), path("${sample_id}*pileup.tsv"), emit: pileup
	tuple val(sample_id), path("${sample_id}*ambiguous.tsv"), emit: ambi
	// path("${sample_id}*stats.tsv"), emit: metrics

	"""
	make_pileup.py ${reference} ${bam} > ${sample_id}.pileup.tsv
	ambiguous_positions.py ${sample_id}.pileup.tsv --output ${sample_id}.ambiguous.tsv
	"""
}