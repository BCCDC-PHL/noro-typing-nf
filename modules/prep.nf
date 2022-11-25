process fastp {

    tag { sample_id }

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
      -t ${task.cpus} 
      -i ${fastq1} \
      -I ${fastq2} \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      -j ${sample_id}.fastp.json \
      --detect_adapter_for_pe
    """
}