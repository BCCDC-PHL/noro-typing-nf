process get_coverage { 
	errorStrategy 'ignore'

	tag { sample_id }

	publishDir "${params.outpath}/${sample_id}/${task.ext.workflow}", pattern: " ${sample_id}*bed", mode:'copy'

	input:
	tuple val(sample_id), path(bamfile), path(bam_index)

	output:
	tuple val(sample_id), path("*bed"), emit: main
	tuple val(sample_id), path("${sample_id}_samtools_provenance.yml"), emit: provenance
	tuple val(sample_id), env(COVERAGE), emit: metric

	"""
	printf -- "- process_name: samtools\\n" > ${sample_id}_samtools_provenance.yml
	printf -- "  tool_name: samtools\\n tool_version: \$(samtools | sed -n '1 p'))\\n" >> ${sample_id}_samtools_provenance.yml

	samtools depth -a ${bamfile} > ${sample_id}_coverage.bed
	COVERAGE=`cat ${sample_id}_coverage.bed | awk '{if(\$3 > ${params.consensus_min_depth} ){total += 1}} END {print total*100/NR}'`
	"""
}

process plot_coverage { 
	errorStrategy 'ignore'

	tag { sample_id }

	conda "${projectDir}/environments/plot.yaml"

	publishDir "${params.outpath}/${sample_id}/${task.ext.workflow}", mode:'copy'

	input:
	tuple val(sample_id), path(coverage_bed)

	output:
	tuple val(sample_id), path("${sample_id}*.pdf")

	"""
	plotting.R ${coverage_bed} ${sample_id}_${task.ext.workflow}_coverage.pdf
	"""
}

process make_pileup {
	errorStrategy 'ignore'

	tag { sample_id }

	publishDir "${params.outpath}/${sample_id}/${task.ext.workflow}/pileups/", pattern: "*{ambiguous.tsv,pileup.tsv}", mode:'copy'

	input:
	tuple val(sample_id), path(reference), path(bam), path(bam_index)

	output:
	tuple val(sample_id), path("${sample_id}*pileup.tsv"), emit: pileup
	tuple val(sample_id), path("${sample_id}*ambiguous.tsv"), emit: ambi
	// path("${sample_id}*stats.tsv"), emit: metrics

	script: 
	"""
	make_pileup.py ${reference} ${bam} > ${sample_id}.pileup.tsv
	ambiguous_positions.py ${sample_id}.pileup.tsv --output ${sample_id}.ambiguous.tsv
	"""
}