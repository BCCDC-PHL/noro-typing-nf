process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/fastp", pattern: "${sample_id}*_R{1,2}.trim.fastq.gz" , mode:'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
    path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastp\\n" > ${sample_id}_fastp_provenance.yml
    printf -- "  tool_name: fastp\\n  tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml
    fastp \
      -t ${task.cpus} \
      -i ${fastq1} \
      -I ${fastq2} \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      -j ${sample_id}.fastp.json \
      --detect_adapter_for_pe
    """
}

process cutadapt {

    tag { sample_id }

    publishDir "${params.outdir}/cutadapt", pattern: "${sample_id}_R{1,2}.trimmed.fastq.gz|${sample_id}.cutadapt.log", mode:'copy'


    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(primers)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: primer_trimmed_reads
    path("${sample_id}.cutadapt.log"), emit: log
    tuple val(sample_id), path("${sample_id}_cutadapt_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: cutadapt\\n" > ${sample_id}_cutadapt_provenance.yml
    printf -- "  tool_name: cutadapt\\n  tool_version: \$(cutadapt --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_cutadapt_provenance.yml
    cutadapt \
      -j ${task.cpus} \
      -g file:${primers} \
      -G file:${primers} \
      -o ${sample_id}_R1.trimmed.fastq.gz \
      -p ${sample_id}_R2.trimmed.fastq.gz \
      ${reads_1} \
      ${reads_2} \
      > ${sample_id}.cutadapt.log
    """
}

process fastQC {

    tag { sample_id }

    publishDir path: "${params.outdir}/fastQC/data", pattern: "${sample_id}_R{1,2}*fastqc.zip", mode: "copy"
    publishDir path: "${params.outdir}/fastQC/", pattern: "${sample_id}_R{1,2}*fastqc.html", mode: "copy"

    input:
    tuple val(sample_id), path(forward), path(reverse)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}*_fastqc.html"), emit: html
    tuple val(sample_id), path("${sample_id}_R{1,2}*_fastqc.zip"), emit: zip
    tuple val(sample_id), path("${sample_id}_fastqc_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastqc\\n" > ${sample_id}_fastqc_provenance.yml
    printf -- "  tool_name: fastqc\\n  tool_version: \$(fastqc --version 2>&1 | sed -n '1 p')\\n" >> ${sample_id}_fastqc_provenance.yml
    fastqc --threads ${task.cpus} ${sample_id}_R1.trim.fastq.gz ${sample_id}_R2.trim.fastq.gz 
    """
}

process fastq_check {

    tag { sample_id }

    conda "${projectDir}/environments/fastx.yaml"

    //publishDir path: "${params.outdir}/fastq_qual/", pattern: "${sample_id}.fqqual.tsv", mode: "copy"
    publishDir path: "${params.outdir}/fastq_qual/", pattern: "${sample_id}.R{1,2}.raw.tsv", mode: "copy"

    input: 
    tuple val(sample_id), path(forward), path(reverse)


    output:
    path("${sample_id}.R{1,2}.raw.tsv"), emit: raw
    path("${sample_id}.QC_fastq.tsv"), emit: formatted

    """
    qc_fastq.sh $forward $reverse
    """
}


process mapping_check {

    tag { sample_id }
    
    //publishDir path: "${params.outdir}/mapping_qual/", pattern: "${bam.simpleName}.mapqual.tsv",mode: "copy"
    publishDir path: "${params.outdir}/mapping_qual/", pattern: "${bam.simpleName}.raw.tsv", mode: "copy"

    input: 
	tuple val(sample_id), path(bam)
	
	output:
    path("${bam.simpleName}.raw.tsv"), emit: raw
    path("${bam.simpleName}.QC_mapping.tsv"), emit: formatted


    """
    qc_mapping.sh $bam 
    """
}

process run_kraken {

    tag {sample_id}

    label 'heavy'

    conda "${projectDir}/environments/kraken.yaml"

    publishDir path: "${params.outdir}/kraken/reports", pattern: "${sample_id}.kraken.report", mode: "copy"
    publishDir path: "${params.outdir}/kraken/output", pattern: "${sample_id}.kraken.out", mode: "copy"


    input: 
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}.kraken.report"), path("${sample_id}.kraken.out")

    """
    kraken2 --confidence 0.1 \
    --threads ${task.cpus} --db ${params.kraken_db} \
    --paired ${reads_1} ${reads_2} \
    --report ${sample_id}.kraken.report > ${sample_id}.kraken.out 
    """
}

process kraken_filter {

    tag {sample_id}

    conda "${projectDir}/environments/kraken.yaml"

    publishDir path: "${params.outdir}/kraken/filtered", pattern: "${sample_id}*fastq.gz", mode: "copy"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(kraken_report), path(kraken_out)

    output:
    tuple val(sample_id), path("${sample_id}_R1.kfilter.fastq.gz"), path("${sample_id}_R2.kfilter.fastq.gz"), emit: filtered

    script:
    """
    # Extract the norovirus specific reads 
    # will exclude anything not under Norovirus clade; including human reads 
    extract_kraken_reads.py \
    -k ${sample_id}.kraken.out -r ${sample_id}.kraken.report \
    -1 ${reads_1} -2 ${reads_2} \
    -o ${sample_id}_R1.kfilter.fastq -o2 ${sample_id}_R2.kfilter.fastq \
    -t ${params.krk_norovirus_id} --include-children &&
    gzip -c ${sample_id}_R1.kfilter.fastq > ${sample_id}_R1.kfilter.fastq.gz &&
    gzip -c ${sample_id}_R2.kfilter.fastq > ${sample_id}_R2.kfilter.fastq.gz

    # # for extracting non-norovirus reads
    # --exclude
    """

}

process run_centrifuge {
    tag {sample_id}

    label 'heavy'

    conda "${projectDir}/environments/kraken.yaml"

    publishDir path: "${params.outdir}/centrifuge/reports", pattern: "${sample_id}.kraken.report", mode: "copy"
    publishDir path: "${params.outdir}/centrifuge/output", pattern: "${sample_id}.kraken.out", mode: "copy"


    input: 
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}.kraken.report"), path("${sample_id}.kraken.out")

    """
    centrifuge  --threads ${task.cpus} \
    -X ${params.centrifuge_db} \
    -1 ${reads_1} -2 ${reads_2} \
    --report-file ${sample_id}.cent.report -S ${sample_id}.cent.out 
    """

}


process build_composite_reference {
    storeDir "${projectDir}/cache/composite"

    input:
    tuple path(human_ref), path(virus_ref)

    output:
    path("${composite_ref}"), emit: fasta
    // path("${composite_ref}.*"), emit: index
    // path("reference_headers.txt"), emit: headers

    script:
    composite_ref = params.composite_ref_name 
    """ 
    cat ${human_ref} ${virus_ref} > ${composite_ref} 
    # bwa index ${composite_ref} 
    # grep ">" ${virus_ref} | cut -c2- > reference_headers.txt
    """
}

process index_composite_reference {
    storeDir "${projectDir}/cache/composite"

    // label 'ultra'

    input:
    path(composite_ref)

    output:
    // val("${composite_ref.simpleName}"), emit: name
    // tuple path("${composite_ref}.bwt"), path("${composite_ref}.amb"), path("${composite_ref}.ann"), path("${composite_ref}.pac"), path("${composite_ref}.sa"), emit: files
    tuple path("${composite_ref}"), path("${composite_ref}.*")

    script:
    ref_name = composite_ref.simpleName
    """
    # bowtie2-build --threads ${task.cpus} ${composite_ref} ${ref_name}
    bwa index -a bwtsw -b 10000000 ${composite_ref}
    """
}

process get_reference_headers {

    input:
    path(viral_reference)

    output:
    path("${ref_name}_headers.txt")

    script:
    ref_name = viral_reference.simpleName

    """
    grep ">" ${viral_reference} | cut -c2- > ${ref_name}_headers.txt
    """    
}

process dehost_fastq {

    label 'heavy'

    tag {sample_id}

    publishDir path: "${params.outdir}/dehosted/", pattern: "${sample_id}*fastq.gz", mode: "copy"
    publishDir path: "${params.outdir}/dehosted/bam", pattern: "${sample_id}*bam", mode: "copy"
    publishDir path: "${params.outdir}/dehosted/metrics", pattern: "${sample_id}*_metrics.txt", mode: "copy"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)
    path(virus_reference_list)
    tuple path(index_name), path("*")

    output:
    tuple val(sample_id), path("${sample_id}_R1.dehost.fastq.gz"), path("${sample_id}_R2.dehost.fastq.gz"), emit: fastq
    tuple val(sample_id), path("${sample_id}.dehosted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}_metrics.txt"), emit: metrics

    script: 
    // params.composite_ref_name
    // dehost.py -k ${virus_names.join(' ')}
    """
    bwa mem -t ${task.cpus} -T 30 ${index_name} ${reads_1} ${reads_2} | \
    dehost.py -r ${virus_reference_list} -n ${sample_id} -o ${sample_id}.dehosted.bam 2> ${sample_id}_metrics.txt
    samtools sort -n --threads ${task.cpus} ${sample_id}.dehosted.bam | \
    samtools fastq --threads ${task.cpus} -1 ${sample_id}_R1.dehost.fastq.gz -2 ${sample_id}_R2.dehost.fastq.gz -0 /dev/null -s /dev/null -n
    """


}