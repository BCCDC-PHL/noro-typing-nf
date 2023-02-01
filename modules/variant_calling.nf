process run_mpileup {

	tag {sample_id}

    publishDir "${params.outdir}/vcf/mpileup", pattern: "${sample_id}.mp.vcf.gz" , mode:'copy'

    input: 
    tuple val(sample_id), path(bamfile), path(bam_index), path(reference), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.mp.vcf.gz")


	"""
	bcftools mpileup --threads ${task.cpus} -q ${params.mp_min_map_qual} -Q ${params.mp_min_base_qual} -m ${params.mp_gapped_indel} -Ou -f ${reference} ${bamfile} | bcftools call --ploidy 1 -Mmv -Oz -o ${sample_id}.mp.vcf
    gzip ${sample_id}.mp.vcf
	"""
}

process run_freebayes {

    tag { sample_id }

    publishDir "${params.outdir}/vcf/freebayes", pattern: "${sample_id}.fb.vcf.gz", mode: 'copy'

    input: 
	tuple val(sample_id), path(bamfile), path(bam_index), path(reference), path(index)
	
	output:
    tuple val(sample_id), path("${sample_id}.fb.vcf.gz")


    """
    freebayes -b ${bamfile} -v ${sample_id}.fb.vcf --gvcf -f ${reference} -P ${params.fb_poly_prob} -p ${params.ploidy} -X -O -m ${params.fb_min_map_qual} -q ${params.fb_min_base_qual} --read-snp-limit ${params.fb_read_snp_limit} -x ${params.fb_indel_excl_window} -F ${params.fb_min_alt_frac} --min-coverage ${params.fb_min_coverage} 
    gzip ${sample_id}.fb.vcf
    """

}

process get_common_snps {

    tag { sample_id }

    publishDir "${params.outdir}/vcf/common", pattern: "${sample_id}.common.vcf.gz", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fb_vcf), path(mp_vcf)
	
    output:
    tuple val(sample_id), path("${sample_id}.common.vcf.gz"), path("${sample_id}*common.vcf.gz.csi")
    
    """
    bcftools sort --output-type b ${fb_vcf} > ${fb_vcf}b
    bcftools sort --output-type b ${mp_vcf} > ${mp_vcf}b

    bcftools index ${fb_vcf}b
    bcftools index ${mp_vcf}b
    
    bcftools isec -n=2 -c all -w1 -O z -o ${sample_id}.common.vcf.gz ${fb_vcf}b ${mp_vcf}b 
    bcftools index ${sample_id}.common.vcf.gz
    """
}