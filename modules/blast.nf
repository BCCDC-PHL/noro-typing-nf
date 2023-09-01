#!/usr/bin/env nextflow

// possibly useful template  https://www.nextflow.io/docs/latest/process.html
// using storeDir will only run the process if the output files DO NOT exist under the specified path
// if they exist, process is skipped and output is passed as the existing files
nextflow.enable.dsl=2

import java.nio.file.Paths

process prep_database {
    storeDir "${projectDir}/cache/blast_db/full/${custom_dir}"
    
    input:
    path(database)

    output:
    path("*filter.fasta")

    script:
    workflow = task.ext.workflow ?: 'full'
    custom_dir = task.ext.custom_dir ?: 'global'
    """
    filter_fasta.py main \
    --header_delim ${params.header_delim} \
    --header_pos_accno ${params.header_pos_accno} \
    --header_pos_type ${params.header_pos_type} \
    --rename \
    ${database} ${database.simpleName}_filter.fasta
    """
}

process extract_genes_blast {

    storeDir "${projectDir}/cache/blast_db/gene/${custom_dir}"

    label 'ultra'
    
    input:
    path(blast_db_fasta)

    output:
    // Format A: used for gene-specific workflows 
    path("*gene*.fasta"), emit: db
    path("*gene*.yaml"), emit: pos
    // Format B: used for global / full-length workflow
    path("*vp1.fasta"), emit: vp1, optional: true
    path("*rdrp.fasta"), emit: rdrp, optional: true
    path("*vp1_pos.yaml"), emit: pos_vp1, optional: true
    path("*rdrp_pos.yaml"), emit: pos_rdrp, optional: true

    script:
    custom_dir = task.ext.custom_dir ?: 'global'
    gene = task.ext.gene ?: 'all'
    blastname = "${blast_db_fasta.simpleName}"
    outname = task.ext.gene ? "${blastname}_gene_${gene}" : "${blastname}_gene_{gene}"

    """
    extract_genes.py database \
    --nthreads ${task.cpus} \
    --query ${blast_db_fasta} \
    --ref ${params.g1_reference} \
    --positions ${params.g1_gene_positions} \
    --gene ${gene} \
    --outfasta ${outname}.fasta \
    --outyaml ${outname}_pos.yaml
    """
}

process make_blast_database {
    storeDir "${projectDir}/cache/blast_db/gene/${custom_dir}"
    
    input:
    path(fasta)

    output:
    path("${db_name}*")

    script:
    workflow = task.ext.workflow ?: ''
    custom_dir = task.ext.custom_dir ?: 'global'
    db_name = task.ext.workflow ? "${workflow}_blastdb.fasta" : "global_database.fasta"

    """
    printf -- "- process_name: blast\\n" > blast_provenance.yml
    printf -- "  tool_name: blastn\\n  tool_version: \$(blastn -version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> blast_provenance.yml

    makeblastdb -dbtype nucl -in ${fasta} -out ${db_name}
    cp ${fasta} ${db_name}
    """
}

process run_self_blast {
    storeDir "${projectDir}/cache/blast_db/gene/${custom_dir}"
    
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
    publishDir "${params.outdir}/blastn/${workflow}/final_refs", pattern: "${sample_id}*fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(contig_file), path(sequence_source)          // 1) Name  2) Assembled contigs  3) FASTA file where the final output sequences are stored
    path(blast_db)                                                          // BLAST database used to search for best matching reference 
    path(self_blast_scores)

    output:
    tuple val(sample_id), path("${sample_id}*blastn.tsv"), emit: raw
    tuple val(sample_id), path("${sample_id}*filter.tsv"), path("${sample_id}*fasta"), emit: main
    //tuple val(sample_id), path("${sample_id}*full.tsv"), emit: full
    tuple val(sample_id), path("${sample_id}*_ref.fasta"), emit: ref
    path("${sample_id}*filter.tsv"), emit: filter

    script:
    contig_mode = params.assemble ? "--contig_mode" : "" 
    workflow = task.ext.workflow ?: 'global'
    blast_db_name = blast_db[0]
    self_blast = self_blast_scores.name != 'NO_FILE' ? "--ref_scores ${self_blast_scores}" : ''
    blast_metric = workflow == 'global' ? params.blast_refsearch_metric : params.blast_typing_metric
    
    """
    blastn -num_threads ${task.cpus} -query ${contig_file} -db ${blast_db_name} -outfmt "6 ${params.blast_outfmt}" > ${sample_id}_${workflow}_blastn.tsv &&
    filter_alignments.py ${sample_id}_${workflow}_blastn.tsv \
    --metric ${blast_metric} \
    ${contig_mode} \
    ${self_blast} \
    --seqs ${sequence_source} \
    --min_id ${params.min_blast_id} \
    --tsv_out ${sample_id}_${workflow}_blastn_filter.tsv \
    --fasta_out ${sample_id}_${workflow}_ref.fasta 

    if [ ${workflow} == 'global' ] ; then
        delim=${params.header_delim}
        pos_type=${params.header_pos_type}
        pos_accno=${params.header_pos_accno}

        header=`head -n1 ${sample_id}_${workflow}_ref.fasta | cut -c2-`
        fields=(\${header//\$delim/ })

        newheader=">${sample_id}|\${fields[\$pos_type]}|\${fields[\$pos_accno]}|sample|global"
        echo \$newheader

        sed -i 1d ${sample_id}_${workflow}_ref.fasta 
        sed -i 1i"\$newheader" ${sample_id}_${workflow}_ref.fasta 
    fi
    """
}

process combine_references {
    
    tag {sample_id}

    errorStrategy 'ignore'

    publishDir "${params.outdir}/blastn/final/", pattern: "${sample_id}*fasta" , mode:'copy'

    input:
    tuple val(sample_id), path(gtype_blast), path(gtype_ref), path(ptype_blast), path(ptype_ref)

    output:
    tuple val(sample_id), path("${sample_id}.ref.final.fasta")
    // path("${sample_id}*final.tsv"), emit: blast

    script: 
    """
    delim=${params.header_delim}
    pos_type=${params.header_pos_type}
    pos_accno=${params.header_pos_accno}

    gheader=`head -n1 ${gtype_ref} | cut -c2-`
    pheader=`head -n1 ${ptype_ref} | cut -c2-`

    if [ -z \$gheader ] && [ -z \$pheader ]; then 
        echo "ERROR: Neither gtype or ptype reference is a valid FASTA"
        exit 1  
    fi

    function concat_header {
        filename=\$1
        append=\$2
        newfile=\$3

        old_header=`head -n1 \$filename`
        full_header="\$old_header \$append"

        echo \$full_header > \$newfile
        sed 1d \$filename >> \$newfile
    }

    gfields=(\${gheader//\$delim/ })
    pfields=(\${pheader//\$delim/ })

    newheader="${sample_id}|\${gfields[\$pos_type]}_\${pfields[\$pos_type]}|\${gfields[\$pos_accno]}_\${pfields[\$pos_accno]}|sample|composite"
    echo \$newheader

    if [ ! -z \$gheader ]; then 
        concat_header ${gtype_ref} \${newheader} gtype_rename.fasta
    fi
    if [ ! -z \$pheader ]; then 
        concat_header ${ptype_ref} \${newheader} ptype_rename.fasta
    fi

    if [ -z \$gheader ] && [ ! -z \$pheader ]; then
        cp ptype_rename.fasta ${sample_id}.ref.final.fasta
    elif [ -z \$pheader ] && [ ! -z \$gheader ]; then
        cp gtype_rename.fasta ${sample_id}.ref.final.fasta
    elif [ \${gfields[\$pos_accno]} == \${pfields[\$pos_accno]} ] ; then 
        cp gtype_rename.fasta ${sample_id}.ref.final.fasta
    else
        cat gtype_rename.fasta ptype_rename.fasta > ${sample_id}.ref.final.fasta
    fi 
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