#!/usr/bin/env nextflow

// possibly useful template  https://www.nextflow.io/docs/latest/process.html
// using storeDir will only run the process if the output files DO NOT exist under the specified path
// if they exist, process is skipped and output is passed as the existing files
nextflow.enable.dsl=2

import java.nio.file.Paths

process make_blast_database {
    storeDir "${projectDir}/cache/blast_db/"
    
    input:
    path(fasta)

    output:
    tuple path("${fasta}"), path("${fasta}.*")

    script:
    """
    makeblastdb -dbtype nucl -in ${fasta} -out ${fasta}
    """
}

process run_blastn {

    tag {sample_id}

    // publishDir "${params.outdir}/blastn", pattern: "${sample_id}_blastn.tsv" , mode:'copy'

    input:
    tuple val(sample_id), path(contig_file)
    tuple path(blast_db), path("*")

    output:
    tuple val(sample_id), path("${sample_id}_blastn.tsv")

    script:
    """
    blastn -query ${contig_file} -db ${blast_db} -outfmt "6 qseqid sseqid pident bitscore qlen slen" > ${sample_id}_blastn.tsv
    """
}

process filter_alignments {

    tag {sample_id}

    publishDir "${params.outdir}/blastn/", pattern: "${sample_id}_blastn_filter.tsv" , mode:'copy'

    input: 
    tuple val(sample_id), path(blast_output)

    output:
    tuple val(sample_id), path("${sample_id}_blastn_filter.tsv")

    script:
    """
    filter_alignments.py ${blast_output} --min_cov ${params.min_blast_cov} --min_id ${params.min_blast_id} --output "${sample_id}_blastn_filter.tsv"
    """

}

process get_best_references {

    tag {sample_id}

    publishDir "${params.outdir}/blastn/final_refs", pattern: "${sample_id}.ref.fasta" , mode:'copy'

    input:
    tuple val(sample_id), path(blast_filtered), path(references)

    output:
    tuple val(sample_id), path("${sample_id}.ref.fasta")

    script: 
    contig_mode = params.assemble ? "-c" : "" 

    """
    get_references.py -b ${blast_filtered} -s ${references} ${contig_mode} -o ${sample_id}.ref.fasta
    """
}