process assembly {

	errorStrategy 'ignore'

	tag {sample_id}

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}_contigs.fa", mode:'copy'
	//publishDir "${custom_outdir}", pattern: "${sample_id}.spades.tar.gz"

	input:
	tuple val(sample_id), path(reads_1), path(reads_2)

	output: 
	tuple val(sample_id), path("${sample_id}_contigs.fa"), emit: contigs
	tuple val(sample_id), path("${sample_id}_assembly_provenance.yml"), emit: provenance

	script:
	
	"""
	printf -- "- process_name: assembly\\n" > ${sample_id}_assembly_provenance.yml
    printf -- "  tool_name: spades\\n  tool_version: \$(rnaviralspades.py --version 2>&1 | cut -d' ' -f4)\\n" >> ${sample_id}_assembly_provenance.yml

	rnaviralspades.py -1 ${reads_1} -2 ${reads_2} -o ${sample_id} &&
	cp ${sample_id}/contigs.fasta ./${sample_id}_contigs.fa &&
	tar -czvf ${sample_id}.spades.tar.gz ${sample_id}/*
	"""

}