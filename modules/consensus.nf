process mask_low_coverage {
	tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}", pattern: "${sample_id}*.bed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bamfile), path(bam_index)
	
    output:
    tuple val(sample_id), path("${output_name}_low_coverage.bed"), emit: bed
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance


    script:
    output_name = "${sample_id}_${task.ext.workflow}"
    """
    printf -- "- process_name: mask_low_coverage\\n" > ${sample_id}-bedtools-provenance.yml
    printf -- "  tools: \\n  - tool_name: bedtools\\n    tool_version: \$(bedtools --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}-bedtools-provenance.yml

    bedtools genomecov -bga -ibam ${bamfile} |  awk '\$4 < ${params.consensus_min_depth} {{print}}' | awk 'BEGIN{FS=OFS="\\t"} {print \$1,\$2+1,\$3+1,\$4}' > ${output_name}_low_coverage.bed
    """
}

process make_consensus {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}", pattern: "${sample_id}*consensus.fasta", mode: 'copy'
    
    input:
    tuple val(sample_id), path(common_vcf), path(vcf_index), path(reference), path(mask_file)
	
    output:
    tuple val(sample_id), path("${output_name}*consensus.fasta"), emit: consensus
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
    output_name = "${sample_id}_${task.ext.workflow}"
    """
    printf -- "- process_name: make_consensus\\n" > ${sample_id}-bcftools-provenance.yml
    printf -- "  tools: \\n  - tool_name: bcftools\\n    tool_version: \$(bcftools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-bcftools-provenance.yml

    bcftools consensus -m ${mask_file} -f ${reference} ${common_vcf} > ${output_name}_consensus.fasta
    """
}