process multiqc {

    conda 'multiqc'

    publishDir "${params.outdir}/multi_qc", pattern: "*report.html", mode:'copy'

    input:
    path('*')

    output:
    path('*report.html')

    script:
    """
    multiqc . -n ${params.run_name}_multiqc_report.html
    
    """
}
