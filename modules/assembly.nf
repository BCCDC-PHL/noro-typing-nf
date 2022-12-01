process assembly {

	tag {sample_id}

    publishDir "${params.outdir}/assembly", pattern: "${sample_id}_contigs.fa", mode:'copy'

	input:
	tuple val(sample_id), path(reads_1), path(reads_2)

	output: 
	tuple val(sample_id), path("${sample_id}_contigs.fa")
	
	"""
	rnaviralspades.py -1 ${reads_1} -2 ${reads_2} -o ${sample_id}_contigs.fa
	"""

}