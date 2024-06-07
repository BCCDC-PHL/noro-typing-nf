#!/usr/bin/env nextflow

// possibly useful template  https://www.nextflow.io/docs/latest/process.html
// using storeDir will ONLY run the process if the output files do not exist under the specified path
// if they exist, process is skipped and output is passed as the existing files
nextflow.enable.dsl=2

import java.nio.file.Paths

process prep_database {

    storeDir "${params.pipeline_cache}/${custom_dir}"
    
    input:
    path(database)

    output:
    path("*filter.fasta")

    script:

    custom_dir = [task.ext.workflow, task.ext.subworkflow].join("/")

    """
    filter_fasta.py \
    --header_delim ${params.header_delim} \
    --header_pos_accno ${params.header_pos_accno} \
    --header_pos_type ${params.header_pos_type} \
    --rename \
    ${database} ${database.simpleName}_filter.fasta
    """
}

process extract_genes_database {

    storeDir "${params.pipeline_cache}/${custom_dir}"

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
    gene = task.ext.gene ?: 'all'
    custom_dir = [task.ext.workflow, task.ext.subworkflow].join("/")
    output_name = task.ext.gene != 'all' ? "${blast_db_fasta.simpleName}_gene_${gene}" : "${blast_db_fasta.simpleName}_gene_{gene}"

    """
    extract_genes.py \
    --nthreads ${task.cpus} \
    --query ${blast_db_fasta} \
    --ref ${params.g1_reference} \
    --positions ${params.g1_gene_positions} \
    --gene ${gene} \
    --outfasta ${output_name}.fasta \
    --outyaml ${output_name}_pos.yaml
    """
}

process make_blast_database {
    storeDir "${params.pipeline_cache}/${custom_dir}"
    
    input:
    path(fasta)

    output:
    path("${blast_db_name}*")

    script:
    custom_dir = [task.ext.workflow, task.ext.subworkflow].join("/")
    blast_db_name = "${task.ext.workflow}_${task.ext.subworkflow}_blastdb.fasta"
    
    """
    printf -- "- process_name: blast\\n" > blast_provenance.yml
    printf -- "  tool_name: blastn\\n  tool_version: \$(blastn -version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> blast_provenance.yml

    makeblastdb -dbtype nucl -in ${fasta} -out ${blast_db_name}
    cp ${fasta} ${blast_db_name}
    """
}

process run_self_blast {
    storeDir "${params.pipeline_cache}/${custom_dir}"
    
    input:
    path(blast_db)

    output:
    path("${outfile}")

    script:
    blast_db_name = blast_db[0]
    custom_dir = [task.ext.workflow, task.ext.subworkflow].join("/")
    outfile = "${blast_db_name.simpleName}_self_blast_scores.tsv"
    """
    blastn -db ${blast_db_name} -query ${blast_db_name} -outfmt "6 qseqid sseqid score" > self_blast.tsv
    awk '{ if (\$1 == \$2) print \$1"\t"\$3}' self_blast.tsv > ${outfile}
    sed -i 1i"name\tselfscore" ${outfile}
    """
}

process run_blastn {

    errorStrategy 'ignore'

    label 'medium'

    tag {sample_id}

    publishDir "${custom_outdir}", pattern: "${output_name}*blastn.tsv" , mode:'copy'
    publishDir "${custom_outdir}", pattern: "${output_name}*filter.tsv" , mode:'copy'
    publishDir "${custom_outdir}", pattern: "${output_name}*fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(contig_file), path(sequence_source)          // 1) Name  2) Assembled contigs  3) FASTA file where the final output sequences are stored
    path(blast_db)                                                          // BLAST database used to search for best matching reference 
    path(self_blast_scores)

    output:
    tuple val(sample_id), path("${output_name}*blastn.tsv"), emit: raw
    tuple val(sample_id), path("${output_name}*filter.tsv"), path("${sample_id}*fasta"), emit: main
    tuple val(sample_id), path("${output_name}*_ref.fasta"), emit: ref
    path("${sample_id}*filter.tsv"), emit: filter

    script:
    workflow = task.ext.workflow
    contig_mode = params.assemble ? "--contig_mode" : "" 
    blast_db_name = blast_db[0]
    self_blast = self_blast_scores.name != 'NO_FILE' ? "--self_blast_scores ${self_blast_scores}" : ''
    blast_metric = workflow == 'global' ? params.blast_metrics_global : params.blast_metrics_composite
    output_name = "${sample_id}_${workflow}_${task.ext.subworkflow}"
    custom_outdir = "${params.outpath}/${sample_id}/${workflow}/"
    
    """
    export BLAST_OUTFMT="${params.blast_outfmt}"

    blastn -num_threads ${task.cpus} -query ${contig_file} -db ${blast_db_name} -outfmt "6 ${params.blast_outfmt}" > ${output_name}_blastn.tsv &&

    filter_alignments.py ${output_name}_blastn.tsv \
    --metric ${blast_metric} \
    ${contig_mode} \
    ${self_blast} \
    --seqs ${sequence_source} \
    --min_id ${params.min_blast_id} \
    --tsv_out ${output_name}_blastn_filter.tsv \
    --fasta_out ${output_name}_ref.fasta 

    # need to reformat the global headers here to prep for downstream
    if [ ${workflow} == 'global' ] ; then
        DELIM=${params.header_delim}
        POS_TYPE=${params.header_pos_type}
        POS_ACCNO=${params.header_pos_accno}

        header=`head -n1 ${output_name}_ref.fasta | cut -c2-`
        fields=(\${header//\$DELIM/ })

        NEWHEADER=">${sample_id}|\${fields[\$POS_TYPE]}|\${fields[\$POS_ACCNO]}|global"
        echo \$NEWHEADER

        sed -i 1d ${output_name}_ref.fasta 
        sed -i 1i"\$NEWHEADER" ${output_name}_ref.fasta 
    fi
    """
}

process combine_references {
    
    tag {sample_id}

    errorStrategy 'ignore'

    publishDir "${params.outpath}/${sample_id}/${task.ext.workflow}/", pattern: "${sample_id}*fasta" , mode:'copy'

    input:
    tuple val(sample_id), path(gtype_blast), path(gtype_reference), path(ptype_blast), path(ptype_reference)

    output:
    tuple val(sample_id), path("${output_name}_ref.fasta")
    // path("${sample_id}*final.tsv"), emit: blast

    script: 
    output_name = "${sample_id}_${task.ext.workflow}"
    """
    DELIM=${params.header_delim}
    POS_TYPE=${params.header_pos_type}
    POS_ACCNO=${params.header_pos_accno}

    GHEADER=`head -n1 ${gtype_reference} | cut -c2-`
    PHEADER=`head -n1 ${ptype_reference} | cut -c2-`

    if [ -z \$GHEADER ] && [ -z \$PHEADER ]; then 
        echo "ERROR: Neither gtype or ptype reference is a valid FASTA"
        exit 1  
    fi

    function concat_header {
        FILENAME=\$1
        APPEND=\$2
        OUTFILE=\$3

        old_header=`head -n1 \$FILENAME`
        full_header="\$old_header \$APPEND"

        echo \$full_header > \$OUTFILE
        sed 1d \$FILENAME >> \$OUTFILE
    }

    GFIELDS=(\${GHEADER//\$DELIM/ })
    PFIELDS=(\${PHEADER//\$DELIM/ })

    NEWHEADER="${sample_id}|\${GFIELDS[\$POS_TYPE]}-\${PFIELDS[\$POS_TYPE]}|\${GFIELDS[\$POS_ACCNO]}-\${PFIELDS[\$POS_ACCNO]}|composite"
    echo \$NEWHEADER

    if [ ! -z \$GHEADER ]; then 
        concat_header ${gtype_reference} \${NEWHEADER} gtype_rename.fasta
    fi
    if [ ! -z \$PHEADER ]; then 
        concat_header ${ptype_reference} \${NEWHEADER} ptype_rename.fasta
    fi

    OUTFILE="${output_name}_ref.fasta"

    if [ -z \$GHEADER ] && [ ! -z \$PHEADER ]; then
        cp ptype_rename.fasta \$OUTFILE
    elif [ -z \$PHEADER ] && [ ! -z \$GHEADER ]; then
        cp gtype_rename.fasta \$OUTFILE
    elif [ \${GFIELDS[\$POS_ACCNO]} == \${PFIELDS[\$POS_ACCNO]} ] ; then 
        cp gtype_rename.fasta \$OUTFILE
    else
        cat gtype_rename.fasta ptype_rename.fasta > \$OUTFILE
    fi 
    """
}

process run_blastx {

    tag {sample_id}

    publishDir "${params.outpath}/blastx/${workflow_type}", pattern: "${sample_id}*.tsv" , mode:'copy'

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