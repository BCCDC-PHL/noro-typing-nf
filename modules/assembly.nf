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