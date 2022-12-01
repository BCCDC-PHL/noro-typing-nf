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

    conda "${baseDir}/environments/fastx.yaml"

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

process kraken {

    tag {sample_id}

    conda "${launchDir}/environment/kraken.yaml"

    publishDir path: "${params.outdir}/kraken/", pattern: "${sample_id}.kraken.report", mode: "copy"
    publishDir path: "${params.outdir}/kraken/raw", pattern: "${sample_id}.kraken.raw", mode: "copy"


    input: 
    tuple val(sample_id), path(forward), path(reverse)

    output:
    path("${sample_id}.kraken.report"), emit: report
    path("${sample_id}.kraken.raw"), emit: raw

    """
    kraken2 --threads ${task.cpus} --db ${params.kraken_db} --paired ${forward} ${reverse} --report ${sampleName}.kraken.report > ${sampleName}.kraken.raw
    """
}
