process fastqc {

    tag { sample_id }

    //publishDir path: "${params.outpath}/qc/fastqc_data", pattern: "${sample_id}_R{1,2}*fastqc.zip", mode: "copy"
    publishDir path: "${params.outpath}/qc/fastqc_html", pattern: "${sample_id}_R{1,2}*fastqc.html", mode: "copy"

    input:
    tuple val(sample_id), path(forward), path(reverse)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}*_fastqc.html"), emit: html
    tuple val(sample_id), path("${sample_id}_R1*fastqc.zip"), path("${sample_id}_R2*fastqc.zip"), emit: zip
    tuple val(sample_id), path("${sample_id}_fastqc_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastqc\\n" > ${sample_id}_fastqc_provenance.yml
    printf -- "  tool_name: fastqc\\n  tool_version: \$(fastqc --version 2>&1 | sed -n '1 p')\\n" >> ${sample_id}_fastqc_provenance.yml

    fastqc --threads ${task.cpus} ${sample_id}_R1.trim.fastq.gz ${sample_id}_R2.trim.fastq.gz 
    """
}


process run_qualimap {
    
    tag { sample_id }

    errorStrategy 'ignore'

    conda 'qualimap'
    
    publishDir path: "${params.outpath}/${sample_id}/${task.ext.workflow}/", pattern: "${sample_id}*pdf", mode: "copy"
    //publishDir path: "${params.outpath}/${sample_id}/qc/mapping/", pattern: "${sample_id}_qmap", mode: "copy"

    input:
    tuple val(sample_id), path(bam_file), path(bam_index)

    output:
    //path("${sample_id}_qmap"), emit: main
    path("${sample_id}*pdf"), emit: pdf
    tuple val(sample_id), path("${sample_id}*provenance.yml"),  emit: provenance


    script:
    """
    printf -- "- process_name: run_qualimap\\n" > ${sample_id}_qualimap_provenance.yml
    printf -- "  tool_name: qualimap\\n  tool_version: \$(bwa 2>&1 |  sed -n '3p' | cut -d' ' -f2)\\n" >> ${sample_id}_qualimap_provenance.yml

    # qualimap bamqc -bam ${bam_file} -outdir ${sample_id}_qmap
    qualimap bamqc -bam ${bam_file} -outfile ${sample_id}_qmap.pdf 
    mv ${sample_id}*_stats/${sample_id}*pdf ./${sample_id}_${task.ext.workflow}_qmap.pdf
    """    
}


process run_quast {

    publishDir path: "${params.outpath}/qc", pattern: "quast_plots/*pdf", mode: "copy"
    publishDir path: "${params.outpath}/qc", pattern: "*{pdf,tsv}", mode: "copy"


    input:
    path(contig_files)

    output:
    path("quast_plots/*pdf"), emit: plots
    path("quast_report.tsv"), emit: tsv
    path("quast_report.pdf"), emit: pdf

    script:
    """
    quast -o quast_results ${contig_files}
    mv quast_results/report.tsv quast_report.tsv
    mv quast_results/report.pdf quast_report.pdf
    mv quast_results/basic_stats/ ./quast_plots
    """    
}

process run_custom_qc {
    errorStrategy 'ignore'
    tag { sample_id }

    publishDir "${params.outpath}/${sample_id}/${task.ext.workflow}", pattern: "${sample_id}*.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bam_index), path(ref), path(consensus)

    output:
    path "${sample_id}.qc.csv", emit: csv
    path "${sample_id}*.png", emit: plot

    script:
    """
    qc.py --outfile ${sample_id}.qc.csv --sample ${sample_id} --ref ${ref} --bam ${bam} --consensus ${consensus}
    """
}

process make_typing_report { 
    
    publishDir "${params.outpath}/", pattern: "*tsv", mode: "copy"

    errorStrategy 'ignore'

    input:
    path(sample_list)
    path(custom_qc_all)
    //path(quast_report)
    path(gblast_collected)
    path(pblast_collected)
    path(global_blast_collected)


    output:
    path("typing_report.tsv")       , emit: main
    path("typing_report_top3.tsv")  , emit: top3


    script:
    global_blast = params.blastn_global_database ? "--globalblast ${global_blast_collected}" : ''
    """
    typing_report.py \
    --sample-list ${sample_list} \
    --gblast ${gblast_collected} \
    --pblast ${pblast_collected} \
    $global_blast \
    --qc ${custom_qc_all} \
    --outpath typing_report.tsv \

    """


}