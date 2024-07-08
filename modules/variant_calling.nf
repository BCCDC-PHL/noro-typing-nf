process mpileup {

	tag {sample_id}

	publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}/vcf", pattern: "${output_name}.mp.vcf.gz" , mode:'copy'

    input: 
    tuple val(sample_id), path(bamfile), path(bam_index), path(reference)

    output:
    tuple val(sample_id), path("${output_name}.mp.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${output_name}_mpileup_provenance.yml"), emit: provenance

    script:
    output_name = "${sample_id}_${task.ext.workflow}"
	"""
	printf -- "- process_name: mpileup\\n" >> ${output_name}_mpileup_provenance.yml
    printf -- "  tool_name: samtools\\n  tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${output_name}_mpileup_provenance.yml
    printf -- "  tool_name: bcftools\\n  tool_version: \$(bcftools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${output_name}_mpileup_provenance.yml

    samtools faidx ${reference}
	bcftools mpileup \
    --threads ${task.cpus} \
    -q ${params.mp_min_map_qual} \
    -Q ${params.mp_min_base_qual} \
    -m ${params.mp_gapped_indel} \
    -Ou \
    -f ${reference} ${bamfile} | bcftools call --ploidy 1 -Mmv -Oz -o ${output_name}.mp.vcf

    gzip ${output_name}.mp.vcf
	"""
}

process freebayes {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}/vcf", pattern: "${output_name}.fb.vcf.gz", mode: 'copy'

    input: 
	tuple val(sample_id), path(bamfile), path(bam_index), path(reference)
	
	output:
    tuple val(sample_id), path("${output_name}.fb.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${output_name}_freebayes_provenance.yml"), emit: provenance

    script:
    output_name = "${sample_id}_${task.ext.workflow}"

    """
	printf -- "- process_name: freebayes\\n" >> ${output_name}_freebayes_provenance.yml
    printf -- "  tool_name: samtools\\n  tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${output_name}_freebayes_provenance.yml
    printf -- "  tool_name: freebayes\\n  tool_version: \$(freebayes --version | cut -d' ' -f3)\\n" >> ${output_name}_freebayes_provenance.yml

    samtools faidx ${reference}
    freebayes --gvcf -X -O \
    -b ${bamfile} \
    -v ${output_name}.fb.vcf \
    -f ${reference} \
    -P ${params.fb_poly_prob} \
    -p ${params.ploidy} \
    -m ${params.fb_min_map_qual} \
    -q ${params.fb_min_base_qual} \
    --read-snp-limit ${params.fb_read_snp_limit} \
    -x ${params.fb_indel_excl_window} \
    -F ${params.fb_min_alt_frac} \
    --min-coverage ${params.fb_min_coverage} 

    gzip ${output_name}.fb.vcf
    """

}

process find_common_snps {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}/", pattern: "${output_name}.common.vcf.gz", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fb_vcf), path(mp_vcf)
	
    output:
    tuple val(sample_id), path("${output_name}.common.vcf.gz"), path("${output_name}*common.vcf.gz.csi"), emit: vcf
    tuple val(sample_id), path("${output_name}_find_common_snps_provenance.yml"), emit: provenance
    
    script:
    output_name = "${sample_id}_${task.ext.workflow}"

    """
    printf -- "- process_name: find_common_snps\\n" >> ${output_name}_find_common_snps_provenance.yml
    printf -- "  tool_name: bcftools\\n  tool_version: \$(bcftools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${output_name}_find_common_snps_provenance.yml

    bcftools sort --output-type b ${fb_vcf} > ${fb_vcf}b
    bcftools sort --output-type b ${mp_vcf} > ${mp_vcf}b

    bcftools index ${fb_vcf}b
    bcftools index ${mp_vcf}b
    
    bcftools isec -n=2 -c all -w1 -O z -o ${output_name}.common.vcf.gz ${fb_vcf}b ${mp_vcf}b 
    bcftools index ${output_name}.common.vcf.gz
    """
}