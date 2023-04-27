process make_union_database {
    storeDir "${projectDir}/cache/blast_db/union_db"

    input:
    path(database_1)
    path(database_2)

    output:
    path("*_union_db.fasta"), emit: fasta
    path("${params.virus_name}_union_headers.txt"), emit: headers

    script:
    """
    merge_fastas.py ${database_1} ${database_2} ${params.virus_name}_union_db.fasta
    grep ">" ${params.virus_name}_union_db.fasta | cut -c2- > ${params.virus_name}_union_headers.txt
    """
}

process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/preprocess/fastp", pattern: "${sample_id}*_R{1,2}.trim.fastq.gz" 
    publishDir "${params.outdir}/preprocess/fastp/html", pattern: "${sample_id}*html" , mode:'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html
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
      -h ${sample_id}.fastp.html \
      --detect_adapter_for_pe
    """
}

process fastp_json_to_csv {
    tag { sample_id }

    publishDir path: "${params.outdir}/preprocess/fastp/csv", pattern: "${sample_id}*csv", mode: "copy"

    input: 
    tuple val(sample_id), path(fastp_json)

    output:
    tuple val(sample_id), path("${sample_id}_fastp.csv")

    """
    fastp_json_to_csv.py ${fastp_json} > ${sample_id}_fastp.csv
    """
}

process cutadapt {

    tag { sample_id }

    publishDir "${params.outdir}/preprocess/cutadapt", pattern: "${sample_id}_R{1,2}.trimmed.fastq.gz"
    publishDir "${params.outdir}/preprocess/cutadapt/logs", pattern: "*cutadapt.log", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: trimmed_reads
    path("${sample_id}.cutadapt.log"), emit: log
    tuple val(sample_id), path("${sample_id}_cutadapt_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: cutadapt\\n" > ${sample_id}_cutadapt_provenance.yml
    printf -- "  tool_name: cutadapt\\n  tool_version: \$(cutadapt --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_cutadapt_provenance.yml
    cutadapt \
      -j ${task.cpus} \
      -g file:${params.primers} \
      -G file:${params.primers_rev} \
      -o ${sample_id}_R1.trimmed.fastq.gz \
      -p ${sample_id}_R2.trimmed.fastq.gz \
      ${reads_1} \
      ${reads_2} \
      > ${sample_id}.cutadapt.log
    """
}

process run_kraken {

    tag {sample_id}

    label 'heavy'

    conda "${projectDir}/environments/kraken.yaml"

    publishDir path: "${params.outdir}/filtering/kraken/reports", pattern: "${sample_id}.kraken.report", mode: "copy"
    publishDir path: "${params.outdir}/filtering/kraken/output", pattern: "${sample_id}.kraken.out", mode: "copy"


    input: 
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}.kraken.report"), path("${sample_id}.kraken.out"), emit: main
    path("${sample_id}.kraken.report"), emit: report

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

    publishDir path: "${params.outdir}/filtering/kraken/filtered", pattern: "${sample_id}*fastq.gz"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(kraken_report), path(kraken_out)

    output:
    tuple val(sample_id), path("${sample_id}_R1.kfilter.fastq.gz"), path("${sample_id}_R2.kfilter.fastq.gz"), emit: fastq

    script:
    """
    # Extract the norovirus specific reads 
    # will exclude anything not under Norovirus clade; including human reads 
    extract_kraken_reads.py  --fastq-output \
    -k ${sample_id}.kraken.out -r ${sample_id}.kraken.report \
    -1 ${reads_1} -2 ${reads_2} \
    -o ${sample_id}_R1.kfilter.fastq -o2 ${sample_id}_R2.kfilter.fastq \
    -t ${params.krk_norovirus_id} 0 --include-children &&
    gzip -c ${sample_id}_R1.kfilter.fastq > ${sample_id}_R1.kfilter.fastq.gz &&
    gzip -c ${sample_id}_R2.kfilter.fastq > ${sample_id}_R2.kfilter.fastq.gz
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
    bwa index -a bwtsw -b 10000000 ${composite_ref}
    """
}

process dehost_fastq {

    label 'heavy'

    tag {sample_id}

    publishDir path: "${params.outdir}/dehosted/", pattern: "${sample_id}*fastq.gz"
    publishDir path: "${params.outdir}/dehosted/bam", pattern: "${sample_id}*bam"
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