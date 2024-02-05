#!/usr/bin/env nextflow

// possibly useful template  https://www.nextflow.io/docs/latest/process.html
// using storeDir will only run the process if the output files DO NOT exist under the specified path
// if they exist, process is skipped and output is passed as the existing files
nextflow.enable.dsl=2

import java.nio.file.Paths

process make_blast_database {
    storeDir "${projectDir}/cache/blast_db/${workflow_type}"
    
    input:
    path(fasta)

    output:
    tuple path("${db_name}"), path("${db_name}.*")

    script:
    db_name = "${fasta.simpleName}"
    workflow_type = "${fasta}" =~ /gtype/ ? "gtype" : "ptype" 
    """
    makeblastdb -dbtype nucl -in ${fasta} -out ${db_name}
    cp -P ${fasta} ${db_name}
    """
}

process run_self_blast {
    storeDir "${projectDir}/cache/blast_db/${workflow_type}"
    
    input:
    tuple path(blast_db), path("*")

    output:
    path("${outfile}")

    script:
    outfile = "${blast_db.simpleName}_ref_scores.tsv"
    workflow_type = "${blast_db}" =~ /gtype/ ? "gtype" : "ptype" 
    """
    blastn -db ${blast_db} -query ${blast_db} -outfmt "6 qseqid sseqid score" > self_blast.tsv
    awk '{ if (\$1 == \$2) print \$1"\t"\$3}' self_blast.tsv > out.tmp
    sed 1i"name\trefscore" out.tmp > ${outfile}
    rm out.tmp
    """
}

process run_blastn {

    tag {sample_id}

    publishDir "${params.outdir}/blastn/${workflow_type}/raw", pattern: "${sample_id}*.tsv" , mode:'copy'

    input:
    tuple val(sample_id), path(contig_file)
    tuple path(blast_db), path("*")

    output:
    tuple val(sample_id), path("${sample_id}*.tsv")

    script:
    workflow_type = "${blast_db}" =~ /gtype/ ? "gtype" : "ptype" 
    """
    blastn -query ${contig_file} -db ${blast_db} -outfmt "6 qseqid saccver qlen qstart qend slen sstart send length pident qcovhsp mismatch gaps evalue bitscore staxids sscinames" > ${sample_id}_blastn.tsv
    """
}

process filter_alignments {

    tag {sample_id}

    publishDir "${params.outdir}/blastn/${workflow_type}", pattern: "${sample_id}*filter.tsv" , mode:'copy'

    input: 
    tuple val(sample_id), path(blast_output), path(self_blast_scores)

    output:
    tuple val(sample_id), path("${sample_id}*filter.tsv")

    script:
    workflow_type = "${self_blast_scores}" =~ /gtype/ ? "gtype" : "ptype" 
    """
    filter_alignments.py ${blast_output} --metric bsr --ref_scores ${self_blast_scores} --min_cov ${params.min_blast_cov} --min_id ${params.min_blast_id} --output ${sample_id}_blastn_${workflow_type}_filter.tsv
    """

}

process get_best_references {

    tag {sample_id}

    publishDir "${params.outdir}/blastn/${workflow_type}/final_refs", pattern: "${sample_id}.ref.fasta" , mode:'copy'

    input:
    tuple val(sample_id), path(blast_filtered), path(reference_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.ref.fasta")

    script: 
    contig_mode = params.assemble ? "--contig_mode" : "" 
    workflow_type = "${blast_filtered}" =~ /gtype/ ? "gtype" : "ptype" 
    """
    get_references.py --blast ${blast_filtered} --metric bsr --seqs ${reference_fasta} ${contig_mode} --output ${sample_id}.ref.fasta
    """
}

process run_blastx {

    tag {sample_id}

    publishDir "${params.outdir}/blastx/", pattern: "${sample_id}*.tsv" , mode:'copy'

    input:
    tuple val(sample_id), path(contig_file)

    output:
    tuple val(sample_id), path("${sample_id}*.tsv")

    script:
    //workflow_type = "${diamond_db}" =~ /gtype/ ? "gtype" : "ptype" 
    """
    # diamond blastx --threads ${task.cpus} -d  -q ${contig_file} -o ${sample_id}_blastx.tsv --outfmt "6 qseqid sseqid pident qlen slen bitscore score"
    blastx -remote -db nr -query ${contig_file} -outfmt "6 qseqid sseqid pident length saccver scinames scomnames sblastnames sskingdoms stitle staxids evalue bitscore" > ${sample_id}.blastx.tsv
    """
}