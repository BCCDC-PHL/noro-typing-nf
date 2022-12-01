process multiqc {

    publishDir "${params.outdir}/multi_qc", pattern: "*_multiqc_report.html", mode:'copy'

    input:
    path('*')

    output:
    path('*_multiqc_*')

    script:
    """
    multiqc . -n ${params.run_name}_multiqc_report.html
    
    """
}
