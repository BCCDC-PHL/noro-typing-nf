process run_mpileup {

	tag {sample_id}

	publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}/vcf", pattern: "${output_name}.mp.vcf.gz" , mode:'copy'

    input: 
    tuple val(sample_id), path(bamfile), path(bam_index), path(reference)

    output:
    tuple val(sample_id), path("${output_name}.mp.vcf.gz")

    script:
    output_name = "${sample_id}_${task.ext.workflow}"
	"""
    samtools faidx ${reference}
	bcftools mpileup --threads ${task.cpus} -q ${params.mp_min_map_qual} -Q ${params.mp_min_base_qual} -m ${params.mp_gapped_indel} -Ou -f ${reference} ${bamfile} | bcftools call --ploidy 1 -Mmv -Oz -o ${output_name}.mp.vcf
    gzip ${output_name}.mp.vcf
	"""
}

process run_freebayes {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}/vcf", pattern: "${output_name}.fb.vcf.gz", mode: 'copy'

    input: 
	tuple val(sample_id), path(bamfile), path(bam_index), path(reference)
	
	output:
    tuple val(sample_id), path("${output_name}.fb.vcf.gz")

    script:
    output_name = "${sample_id}_${task.ext.workflow}"

    """
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

process get_common_snps {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/${task.ext.workflow}/", pattern: "${output_name}.common.vcf.gz", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fb_vcf), path(mp_vcf)
	
    output:
    tuple val(sample_id), path("${output_name}.common.vcf.gz"), path("${output_name}*common.vcf.gz.csi")
    
    script:
    output_name = "${sample_id}_${task.ext.workflow}"

    """
    bcftools sort --output-type b ${fb_vcf} > ${fb_vcf}b
    bcftools sort --output-type b ${mp_vcf} > ${mp_vcf}b

    bcftools index ${fb_vcf}b
    bcftools index ${mp_vcf}b
    
    bcftools isec -n=2 -c all -w1 -O z -o ${output_name}.common.vcf.gz ${fb_vcf}b ${mp_vcf}b 
    bcftools index ${output_name}.common.vcf.gz
    """
}