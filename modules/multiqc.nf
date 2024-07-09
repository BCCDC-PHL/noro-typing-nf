process multiqc {

    conda 'multiqc'

    publishDir "${params.outdir}/", pattern: "*report.html", mode:'copy'

    input:
    path('*')

    output:
    path('*report.html')
    // path('provenance.yml'), emit: provenance

    script:
    """
    printf -- "- process_name: multiqc\\n" > multiqc_provenance.yml
    printf -- "  tools: \\n  - tool_name: multiqc\\n    tool_version: \$(multiqc --version 2>&1 | cut -d' ' -f3)\\n" >> multiqc_provenance.yml

    multiqc . --config ${params.multiqc_config} --filename ${params.run_name}_multiqc_report.html
    """
}
