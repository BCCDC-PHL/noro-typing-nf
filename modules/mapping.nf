process map_reads {

    tag {sample_id}

    errorStrategy 'ignore'

	label 'heavy'

    input: 
    tuple val(sample_id), path(reads_1), path(reads_2), path(reference)

    output:
    tuple val(sample_id), path("${output_name}.sam"), emit: sam
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit:provenance

    script:
    output_name = "${sample_id}_${task.ext.workflow}"
	"""
    printf -- "- process_name: map_reads\\n" > ${sample_id}-bwa-provenance.yml
    printf -- "  tools: \\n  - tool_name: bwa\\n    tool_version: \$(bwa 2>&1 | head -n3 | tail -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-bwa-provenance.yml

    bwa index ${reference}
	bwa mem -t ${task.cpus} ${reference} ${reads_1} ${reads_2} > ${output_name}.sam
	"""
}
process sort_filter_sam {

	tag {sample_id}

    //publishDir "${params.outdir}/mapped_reads/sorted${workflow}", pattern: "${sample_id}.sorted.bam"

    input: 
    tuple val(sample_id), path(samfile)

    output:
    tuple val(sample_id), path("${output_name}.sort.filter.bam"), path("${sample_id}*bai"), emit: bam
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
    output_name = "${sample_id}_${task.ext.workflow}"
	"""
    printf -- "- process_name: sort_filter_sam\\n" > ${sample_id}-samtools-provenance.yml
    printf -- "  tools: \\n  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-samtools-provenance.yml

	samtools view -f 3 -F 2828 -q 30 -h ${samfile} | samtools sort -o ${output_name}.sort.filter.bam 
    samtools index ${output_name}.sort.filter.bam 
	"""
}

process merge_fasta_bam {

	tag {sample_id}

    publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}/", pattern: "*bam"
    publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}/", pattern: "*fasta" , mode:'copy'

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
