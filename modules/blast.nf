#!/usr/bin/env nextflow

// possibly useful template  https://www.nextflow.io/docs/latest/process.html
// using storeDir will only run the process if the output files DO NOT exist under the specified path
// if they exist, process is skipped and output is passed as the existing files
nextflow.enable.dsl=2

import java.nio.file.Paths

process make_blast_database {
    //storeDir "${baseDir}/blast_db/"
    
    input:
    path(fasta)

    output:
    tuple val("ref_db"), path("ref_db.*")

    script:
    blast_db = fasta.simpleName
    """
    makeblastdb -dbtype nucl -in ${fasta} -out ref_db
    """
}

process run_blastn {

    tag {sample_id}

    publishDir "${params.outdir}/blastn", pattern: "${sample_id}_blastn.tsv" , mode:'copy'

    input:
    tuple val(sample_id), path(contig_file), val(blast_db), path("*")

    output:
    tuple val(sample_id), path("${sample_id}_blastn.tsv")

    script:
    """
    blastn -query ${contig_file} -db ${blast_db} -outfmt "6 qseqid sseqid pident bitscore qlen slen" > ${sample_id}_blastn.tsv
    """
}