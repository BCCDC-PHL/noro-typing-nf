process mask_low_coverage {
	tag { sample_id }

    publishDir "${params.outdir}/consensus/mask", pattern: "${sample_id}*.bed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bamfile), path(bam_index)
	
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
    tuple val(sample_id), path("${sample_id}.consensus.fasta")

    
    """
    bcftools consensus -m ${mask_file} -f ${reference} ${common_vcf} > ${sample_id}.consensus.fasta &&
    HEADER=`head -n 1 ${sample_id}.consensus.fasta | tr '${params.header_delim}>' '|'` &&
    sed -i 1d ${sample_id}.consensus.fasta &&
    sed -i 1i">${sample_id}\${HEADER}" ${sample_id}.consensus.fasta 
    """
    //# TYPE=`head -n 1 ${sample_id}.consensus.fasta | cut -d"${params.header_delim}" -f${params.header_pos_type},${params.header_pos_strain} --output-delimiter=\|` &&
}