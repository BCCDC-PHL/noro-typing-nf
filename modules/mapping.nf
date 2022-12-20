process create_bwa_index {

    // storeDir "${projectDir}/cache/bwa_index/"
	tag {sample_id}
    
    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), val("${fasta}"), path("${fasta}.*")

    script:
    """
    bwa index ${fasta}
    """

}

process create_fasta_index {

    // storeDir "${projectDir}/cache/bwa_index/"
	tag {sample_id}
    
    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${fasta}"), path("${fasta}.fai")

    script:
    """
	samtools faidx ${fasta}
    """

}

process map_reads {

    tag {sample_id}

	label 'heavy'

    //publishDir "${params.outdir}/aligned_reads/sam", pattern: "${sample_id}.sam" , mode:'copy'

    input: 
    tuple val(sample_id), path(reads_1), path(reads_2), val(index_name), path("*")

    output:
    tuple val(sample_id), path("${sample_id}.sam")


	"""
	bwa mem -t ${task.cpus} ${index_name} ${reads_1} ${reads_2} > ${sample_id}.sam
	"""
}
process sort_filter_sam {

	tag {sample_id}

    publishDir "${params.outdir}/mapped_reads/sorted", pattern: "${sample_id}.sorted.bam" , mode:'copy'

    input: 
    tuple val(sample_id), path(samfile)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

	"""
	samtools view -f 3 -F 2828 -q 30 -h ${samfile} | samtools sort -o ${sample_id}.sorted.bam 
	"""
}


process index_bam {

    publishDir "${params.outdir}/mapped_reads/sorted", pattern: "${sample_id}*.bai" , mode:'copy'

	input:
	tuple val(sample_id), path(bamfile)
	output:
	tuple val(sample_id), path("${sample_id}*.bai")
	"""
	samtools index ${bamfile}
	"""
}