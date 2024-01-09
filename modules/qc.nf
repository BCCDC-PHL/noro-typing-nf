process fastQC {

    tag { sample_id }

    publishDir path: "${params.outdir}/qc/fastQC/data", pattern: "${sample_id}_R{1,2}*fastqc.zip", mode: "copy"
    publishDir path: "${params.outdir}/qc/fastQC/", pattern: "${sample_id}_R{1,2}*fastqc.html", mode: "copy"

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

process fastq_check {

    tag { sample_id }

    conda "${projectDir}/environments/fastx.yaml"

    //publishDir path: "${params.outdir}/fastq_qual/", pattern: "${sample_id}.fqqual.tsv", mode: "copy"
    publishDir path: "${params.outdir}/qc/fastx/", pattern: "${sample_id}.R{1,2}.raw.tsv", mode: "copy"

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
    publishDir path: "${params.outdir}/qc/mapping/", pattern: "${bam.simpleName}.raw.tsv", mode: "copy"

    input: 
	tuple val(sample_id), path(bam)
	
	output:
    path("${bam.simpleName}.raw.tsv"), emit: raw
    path("${bam.simpleName}.QC_mapping.tsv"), emit: formatted


    """
    qc_mapping.sh $bam 
    """
}

process run_qualimap {
    
    tag { sample_id }

    errorStrategy 'ignore'

    conda 'qualimap'
    
    publishDir path: "${params.outdir}/qc/mapping/pdf", pattern: "${sample_id}*pdf", mode: "copy"
    publishDir path: "${params.outdir}/qc/mapping/", pattern: "${sample_id}_qmap", mode: "copy"

    input:
    tuple val(sample_id), path(bam_file), path(bam_index)

    output:
    path("${sample_id}_qmap"), emit: main
    path("${sample_id}.pdf"), emit: pdf

    script:
    """
    qualimap bamqc -bam ${bam_file} -outdir ${sample_id}_qmap &&
    qualimap bamqc -bam ${bam_file} -outfile ${sample_id}.pdf && 
    mv ${sample_id}*_stats/${sample_id}.pdf .
    """    
}


process run_quast {

    publishDir path: "${params.outdir}/qc/assembly/plots", pattern: "quast_results/basic_stats/*pdf", mode: "copy"
    publishDir path: "${params.outdir}/qc/assembly/reports", pattern: "quast_results/*{pdf,html,tsv,log}", mode: "copy"


    input:
    path(contig_files)

    output:
    path("quast_results/basic_stats/*pdf"), emit: plots
    path("quast_results/report.tsv"), emit: tsv
    path("quast_results/report.pdf"), emit: pdf

    script:
    """
    quast -o quast_results ${contig_files}
    """    
}

process run_mapping_qc {
    tag { sample_id }

    publishDir "${params.outdir}/qc/full/plots", pattern: "${sample_id}.depth.png", mode: 'copy'
    publishDir "${params.outdir}/qc/full/raw", pattern: "${sample_id}*csv", mode: 'copy'

    input:
    tuple sample_id, path(bam), path(bam_index), path(consensus), path(ref)

    output:
    path "${sample_id}.qc.csv", emit: csv
    path "${sample_id}.depth.png", emit: plot

    script:

    """
    qc.py ${qcSetting} --outfile ${sample_id}.qc.csv --sample ${sample_id} --ref ${ref} --bam ${bam} --consensus ${consensus}
    """
}

process run_custom_qc {
    errorStrategy 'ignore'
    tag { sample_id }

    publishDir "${params.outdir}/qc/custom/csv", pattern: "${sample_id}*csv", mode: 'copy'
    publishDir "${params.outdir}/qc/custom/plots", pattern: "${sample_id}*png", mode: 'copy'

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
    
    publishDir "${params.outdir}", pattern: "*tsv", mode: "copy"

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
    fullblast = params.reference_db_full ? "--globalblast ${global_blast_collected}" : ''
    """
    typing_report.py \
    --sample-list ${sample_list} \
    --gblast ${gblast_collected} \
    --pblast ${pblast_collected} \
    $fullblast \
    --qc ${custom_qc_all} \
    --outpath typing_report.tsv \

    """


}