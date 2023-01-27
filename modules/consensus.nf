process mask_low_coverage {
	tag { sample_id }

    publishDir "${params.outdir}/consensus/mask", pattern: "${sample_id}*.bed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bamfile)
	
    output:
    tuple val(sample_id), path("${sample_id}_low_coverage.bed")
    
    """
    bedtools genomecov -bga -ibam ${bamfile} |  awk '\$4 < ${params.consensus_min_depth} {{print}}' | awk 'BEGIN{FS=OFS="\\t"} {print \$1,\$2+1,\$3+1,\$4}' > "${sample_id}_low_coverage.bed"
    """
}

process make_consensus {

    tag { sample_id }

    publishDir "${params.outdir}/consensus", pattern: "${sample_id}.consensus.fasta", mode: 'copy'
    
    input:
    tuple val(sample_id), path(common_vcf), path(vcf_index), path(reference), path(mask_file)
	
    output:
    path("${sample_id}.consensus.fasta")
    
    """
    bcftools consensus -m ${mask_file} -f ${reference} ${common_vcf} > ${sample_id}.consensus.fasta &&
    TYPE=`head -n 1 ${sample_id}.consensus.fasta | cut -d"|" -f2,3 --output-delimiter=_` &&
    sed -i 1d ${sample_id}.consensus.fasta &&
    sed -i 1i">${sample_id}_\${TYPE}" ${sample_id}.consensus.fasta 
    """
}