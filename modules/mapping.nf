process map_reads {

    tag {sample_id}

    errorStrategy 'ignore'

	label 'heavy'

    input: 
    tuple val(sample_id), path(reads_1), path(reads_2), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

	"""
    printf -- "- process_name: map_reads\\n" > ${sample_id}_bwa_provenance.yml
    printf -- "  tool_name: bwa\\n  tool_version: \$(bwa 2>&1 | head -n3 | tail -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_bwa_provenance.yml

    bwa index ${reference}
	bwa mem -t ${task.cpus} ${reference} ${reads_1} ${reads_2} > ${sample_id}.sam
	"""
}
process sort_filter_index_sam {

	tag {sample_id}

    //publishDir "${params.outpath}/mapped_reads/sorted${workflow}", pattern: "${sample_id}.sorted.bam"

    input: 
    tuple val(sample_id), path(samfile)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}*bai")

    script:
	"""
    printf -- "- process_name: sort_filter_index_sam\\n" > ${sample_id}_samtools_provenance.yml
    printf -- "  tool_name: samtools\\n  tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_samtools_provenance.yml

	samtools view -f 3 -F 2828 -q 30 -h ${samfile} | samtools sort -o ${sample_id}.sorted.bam 
    samtools index ${sample_id}.sorted.bam 
	"""
}

process merge_fasta_bam {

	tag {sample_id}

    //publishDir "${params.outpath}/mapped_reads/merged", pattern: "*bam"
    publishDir "${params.outpath}/${sample_id}/${task.ext.workflow}/", pattern: "*fasta" , mode:'copy'

    input: 
    tuple val(sample_id),  path(references), path(bam_file), path(bam_index)

    output:
    tuple val(sample_id), path("${sample_id}*fasta"), path("${sample_id}*bam"), path("${sample_id}*bai"), emit: main
    tuple val(sample_id), path("${sample_id}*bam"), path("${sample_id}*bai"), emit: bam
    tuple val(sample_id), path("${sample_id}*fasta"), emit: ref

	"""
    merge_fasta_bam.py ${bam_file} ${references} --outbam ${sample_id}_merged.bam --outfasta ${sample_id}_${task.ext.workflow}_merged_ref.fasta
	"""
}
