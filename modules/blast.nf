#!/usr/bin/env nextflow

// possibly useful template  https://www.nextflow.io/docs/latest/process.html
// using storeDir will only run the process if the output files DO NOT exist under the specified path
// if they exist, process is skipped and output is passed as the existing files
nextflow.enable.dsl=2

import java.nio.file.Paths

process prep_database {
    storeDir "${projectDir}/cache/blast_db/full_${workflow}"
    
    input:
    path(database)

    output:
    path("*filter.fasta")

    script:
    workflow = task.ext.workflow ?: 'full_genome'
    """
    filter_fasta.py \
    --header_delim ${params.header_delim} \
    --header_pos_accno ${params.header_pos_accno} \
    --header_pos_type ${params.header_pos_type} \
    ${database} ${database.simpleName}_filter.fasta
    """
}

process extract_genes_blast {
    storeDir "${projectDir}/cache/blast_db/gene_${custom_dir}"
    
    input:
    path(blast_db)

    output:
    path("${blast_db.simpleName}_${gene}.fasta")

    script:
    custom_dir = task.ext.custom_dir ?: 'full_genome'
    gene = task.ext.gene ?: ''
    """
    extract_genes.py database --query ${blast_db} --ref ${params.g1_reference} --positions ${params.g1_gene_positions} --gene ${gene} --output ${blast_db.simpleName}_${gene}.fasta
    """
}

process make_blast_database {
    storeDir "${projectDir}/cache/blast_db/gene_${custom_dir}"
    
    input:
    path(fasta)

    output:
    path("${db_name}*")

    script:
    workflow = task.ext.workflow ?: ''
    db_name = task.ext.workflow ? "${workflow}_blastdb.fasta" : "reference_db.fasta"

    """
    makeblastdb -dbtype nucl -in ${fasta} -out ${db_name}
    cp ${fasta} ${db_name}
    """
}

process run_self_blast {
    storeDir "${projectDir}/cache/blast_db/gene_${custom_dir}"
    
    input:
    path(blast_db)

    output:
    path("${outfile}")

    script:
    db_name = blast_db[0]
    custom_dir = task.ext.custom_dir ?: 'full_genome'
    outfile = "${db_name.simpleName}_ref_scores.tsv"
    workflow = task.ext.workflow ?: ''
    """
    blastn -db ${db_name} -query ${db_name} -outfmt "6 qseqid sseqid score" > self_blast.tsv
    awk '{ if (\$1 == \$2) print \$1"\t"\$3}' self_blast.tsv > ${outfile}
    sed -i 1i"name\trefscore" ${outfile}
    """
}

process run_blastn {

    errorStrategy 'ignore'

    label 'medium'

    tag {sample_id}

    publishDir "${params.outdir}/blastn/${workflow}/raw", pattern: "${sample_id}*blastn.tsv" , mode:'copy'
    publishDir "${params.outdir}/blastn/${workflow}/filtered", pattern: "${sample_id}*filter.tsv" , mode:'copy'
    publishDir "${params.outdir}/blastn/${workflow}/full", pattern: "${sample_id}*full.tsv" , mode:'copy'
    publishDir "${params.outdir}/blastn/${workflow}/final_refs", pattern: "${sample_id}*fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(contig_file), path(reference_fasta), path(self_blast_scores)
    path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}*blastn.tsv"), emit: raw
    tuple val(sample_id), path("${sample_id}*filter.tsv"), path("${sample_id}*fasta"), emit: main
    tuple val(sample_id), path("${sample_id}*full.tsv"), emit: full
    path("${sample_id}*filter.tsv"), emit: filter

    script:
    contig_mode = params.assemble ? "--contig_mode" : "" 
    workflow = task.ext.workflow ?: ''
    blast_db_name = blast_db[0]
    """
    blastn -num_threads ${task.cpus} -query ${contig_file} -db ${blast_db_name} -outfmt "6 ${params.blast_outfmt}" > ${sample_id}_${workflow}_blastn.tsv &&
    filter_alignments.py ${sample_id}_${workflow}_blastn.tsv \
    --metric prop_covered \
    ${contig_mode} \
    --seqs ${reference_fasta} \
    --ref_scores ${self_blast_scores} \
    --min_id ${params.min_blast_id} \
    --tsv_out ${sample_id}_${workflow}_blastn_filter.tsv \
    --fasta_out ${sample_id}_${workflow}_ref.fasta 
    """
}

process find_reference {

    errorStrategy 'ignore'

    label 'medium'

    tag {sample_id}

    publishDir "${params.outdir}/blastn/${workflow}/raw", pattern: "${sample_id}*blastn.tsv" , mode:'copy'
    publishDir "${params.outdir}/blastn/${workflow}/filtered", pattern: "${sample_id}*filter.tsv" , mode:'copy'
    publishDir "${params.outdir}/blastn/${workflow}/full", pattern: "${sample_id}*full.tsv" , mode:'copy'
    publishDir "${params.outdir}/blastn/${workflow}/final_refs", pattern: "${sample_id}*fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(contig_file)
    path(reference_database)

    output:
    tuple val(sample_id), path("${sample_id}*blastn.tsv"), emit: raw
    tuple val(sample_id), path("${sample_id}*filter.tsv"), path("${sample_id}*fasta"), emit: main
    tuple val(sample_id), path("${sample_id}*full.tsv"), emit: full
    tuple val(sample_id), path("${sample_id}*_ref.fasta"), emit: ref
    path("${sample_id}*filter.tsv"), emit: filter

    script:
    contig_mode = params.assemble ? "--contig_mode" : "" 
    workflow = task.ext.workflow ?: ''

    """
    blastn -num_threads ${task.cpus} -query ${contig_file} -db ${reference_database[0]} -outfmt "6 ${params.blast_outfmt}" > ${sample_id}_${workflow}_blastn.tsv &&
    filter_alignments.py ${sample_id}_${workflow}_blastn.tsv \
    --metric prop_covered \
    ${contig_mode} \
    --seqs ${reference_database[0]} \
    --min_id ${params.min_blast_id} \
    --tsv_out ${sample_id}_db_blastn_filter.tsv \
    --fasta_out ${sample_id}_db_ref.fasta 
    """
}



process combine_references {
    
    tag {sample_id}

    errorStrategy 'ignore'

    publishDir "${params.outdir}/blastn/final/", pattern: "${sample_id}*fasta" , mode:'copy'

    input:
    tuple val(sample_id), path(gtype_blast), path(gtype_ref), path(ptype_blast), path(ptype_ref)

    output:
    tuple val(sample_id), path("${sample_id}.ref.final.fasta"), emit: ref
    path("${sample_id}*final.tsv"), emit: blast

    script: 
    """
    combine_references.py \
    -g ${gtype_blast} \
    -p ${ptype_blast} \
    -G ${gtype_ref} \
    -P ${ptype_ref} \
    --outfasta ${sample_id}.ref.final.fasta \
    --outblast ${sample_id}.blast.final.tsv
    """
}

process run_blastx {

    tag {sample_id}

    publishDir "${params.outdir}/blastx/${workflow_type}", pattern: "${sample_id}*.tsv" , mode:'copy'

    input:
    tuple val(sample_id), path(contig_file)
    tuple path(diamond_db), path("*")

    output:
    tuple val(sample_id), path("${sample_id}*.tsv")

    script:
    workflow_type = "${diamond_db}" =~ /gtype/ ? "gtype" : "ptype" 
    """
    diamond blastx --threads ${task.cpus} -d ${diamond_db} -q ${contig_file} -o ${sample_id}_blastx.tsv --outfmt "6 qseqid sseqid pident qlen slen length bitscore score"
    """
}