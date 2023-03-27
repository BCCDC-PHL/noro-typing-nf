process multiqc {

    conda 'multiqc'

    publishDir "${params.outdir}/multi_qc", pattern: "*report.html", mode:'copy'

    input:
    path('*')

    output:
    path('*report.html')

    script:
    """
    multiqc . --config ${params.multiqc_config} --filename ${params.run_name}_multiqc_report.html
    
    """
}
