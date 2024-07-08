process hash_files {

    tag { sample_id }

    input:
    tuple  val(sample_id), path(files_to_hash), val(file_type)

    output:
    tuple  val(sample_id), path("${sample_id}_${file_type}.sha256.csv"), emit: csv
    tuple  val(sample_id), path("${sample_id}_${file_type}_provenance.yml"), emit: provenance

    script:
    """
    shasum -a 256 ${files_to_hash} | tr -s ' ' ',' > ${sample_id}_${file_type}.sha256.csv
    while IFS=',' read -r hash filename; do
        printf -- "- input_filename: \$filename\\n  input_path: \$(realpath \$filename)\\n  sha256: \$hash\\n" >> ${sample_id}_${file_type}_provenance.yml
    done < ${sample_id}_${file_type}.sha256.csv
    """
}

process make_union_database {

    storeDir "${params.pipeline_cache}/dehost"

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

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*{csv,html}", mode: "copy"

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html
    tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance
    tuple val(sample_id), path("${sample_id}_fastp.csv"), emit: csv

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
      # --umi_loc=per_read --umi_len=8 \
      --detect_adapter_for_pe && 
    fastp_json_to_csv.py ${sample_id}.fastp.json > ${sample_id}_fastp.csv
    """
}

process cutadapt {

    tag { sample_id }

    // publishDir "${params.outdir}/preprocess/cutadapt", pattern: "${sample_id}_R{1,2}.trimmed.fastq.gz"
    publishDir "${params.outdir}/${sample_id}", pattern: "*cutadapt.log", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: trimmed_reads
    path("${sample_id}_cutadapt.log"), emit: log
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
      > ${sample_id}_cutadapt.log
    """
}

process run_kraken {

    tag {sample_id}

    label 'heavy'

    conda "${projectDir}/environments/kraken.yaml"

    publishDir path: "${params.outdir}/${sample_id}", pattern: "${sample_id}*report", mode: "copy"
    publishDir path: "${params.outdir}/${sample_id}", pattern: "${sample_id}*out", mode: "copy"

    input: 
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}*report"), path("${sample_id}*out"), emit: main
    tuple val(sample_id), path("${sample_id}_R1.kfilter.fastq.gz"), path("${sample_id}_R2.kfilter.fastq.gz"), emit: fastq
    path("${sample_id}_kraken.report"), emit: report

    """
    printf -- "- process_name: kraken2\\n" > ${sample_id}_kraken2_provenance.yml
    printf -- "  tool_name: kraken2\\n  tool_version: \$(kraken2 --version | head -n1 | cut -d' ' -f3)\\n" >> ${sample_id}_kraken2_provenance.yml

    kraken2 --confidence 0.1 \
    --threads ${task.cpus} --db ${params.kraken_db} \
    --paired ${reads_1} ${reads_2} \
    --report ${sample_id}_kraken.report > ${sample_id}_kraken.out &&

    # Extract the norovirus specific reads 
    extract_kraken_reads.py  --fastq-output \
    -k ${sample_id}_kraken.out -r ${sample_id}_kraken.report \
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

    publishDir path: "${params.outdir}/centrifuge/reports", pattern: "${sample_id}_kraken.report", mode: "copy"
    publishDir path: "${params.outdir}/centrifuge/output", pattern: "${sample_id}_kraken.out", mode: "copy"


    input: 
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_kraken.report"), path("${sample_id}_kraken.out")

    """
    centrifuge  --threads ${task.cpus} \
    -X ${params.centrifuge_db} \
    -1 ${reads_1} -2 ${reads_2} \
    --report-file ${sample_id}.cent.report -S ${sample_id}.cent.out 
    """

}


process build_composite_reference {
    storeDir "${params.pipeline_cache}/dehost"

    input:
    tuple path(human_ref), path(virus_ref)

    output:
    tuple path("composite_reference.fasta"), path("composite_reference.fasta.*")
    // path("${composite_ref}.*"), emit: index
    // path("reference_headers.txt"), emit: headers

    script: 
    """ 
    echo "Building composite reference..." && 
    cat ${human_ref} ${virus_ref} > composite_reference.fasta &&
    echo "Indexing composite reference..." && 
    bwa index -a bwtsw -b 10000000 composite_reference.fasta &&
    echo "Complete." 
    """
}

process dehost_fastq {

    label 'heavy'

    tag {sample_id}

    // publishDir path: "${params.outdir}/preprocess/dehosted/", pattern: "${sample_id}*fastq.gz"
    // publishDir path: "${params.outdir}/preprocess/dehosted/bam", pattern: "${sample_id}*bam"
    publishDir path: "${params.outdir}/${sample_id}/", pattern: "${sample_id}*_metrics.txt", mode: "copy"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)
    path(virus_reference_list)
    tuple path(index_name), path("*")

    output:
    tuple val(sample_id), path("${sample_id}_R1.dehost.fastq.gz"), path("${sample_id}_R2.dehost.fastq.gz"), emit: fastq
    tuple val(sample_id), path("${sample_id}.dehosted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}_dehost_metrics.txt"), emit: metrics

    script: 
    // params.composite_ref_name
    // dehost.py -k ${virus_names.join(' ')}
    """
    bwa mem -t ${task.cpus} -T 30 ${index_name} ${reads_1} ${reads_2} | \
    dehost.py -r ${virus_reference_list} -n ${sample_id} -o ${sample_id}.dehosted.bam 2> ${sample_id}_dehost_metrics.txt
    samtools sort -n --threads ${task.cpus} ${sample_id}.dehosted.bam | \
    samtools fastq --threads ${task.cpus} -1 ${sample_id}_R1.dehost.fastq.gz -2 ${sample_id}_R2.dehost.fastq.gz -0 /dev/null -s /dev/null -n
    """
}