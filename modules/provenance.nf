
process collect_provenance {

  tag { sample_id }

  executor 'local'

  publishDir params.versioned_outdir ? "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output/Provenance_files" : "${params.outdir}/${params.run_name}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output/Provenance_files", pattern: "${sample_id}_*_provenance.yml", mode: 'copy'

  input:
  tuple val(sample_id), path(provenance_files)

  output:
  tuple val(sample_id), file("${sample_id}_*_provenance.yml")

  script:
  """
  cat ${provenance_files} > ${sample_id}_\$(date +%Y%m%d%H%M%S)_provenance.yml
  """
}
process full_provenance {

	tag { sample_id }

	executor 'local'

	publishDir "${params.outdir}/${params.run_name}-provenance.yml", pattern: "provenance.yml", mode: 'copy'

	input:
	path(provenance_files)

	output:
	file("provenance.yml")

	script:
	"""
    printf -- "- process_name: fastp\\n" >> provenance.yml
    printf -- "  tool_name: fastp\\n  tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> provenance.yml

	printf -- "- process_name: cutadapt\\n" >> provenance.yml
    printf -- "  tool_name: cutadapt\\n  tool_version: \$(cutadapt --version 2>&1 | cut -d ' ' -f 2)\\n" >> provenance.yml

	printf -- "- process_name: kraken2\\n" >> provenance.yml
    printf -- "  tool_name: kraken2\\n  tool_version: \$(kraken2 --version | head -n1 | cut -d' ' -f3)\\n" >> provenance.yml

	printf -- "- process_name: multiqc\\n" >> provenance.yml
    printf -- "  tool_name: multiqc\\n  tool_version: \$(multiqc --version 2>&1 | cut -d' ' -f3)\\n" >> provenance.yml

	printf -- "- process_name: spades\\n" >> provenance.yml
    printf -- "  tool_name: spades\\n  tool_version: \$(rnaviralspades.py --version 2>&1 | cut -d' ' -f4)\\n" >> provenance.yml

	printf -- "- process_name: mask_low_coverage\\n" >> provenance.yml
    printf -- "  tool_name: bedtools\\n  tool_version: \$(bedtools --version 2>&1 | cut -d' ' -f2)\\n" >> provenance.yml

	printf -- "- process_name: make_consensus\\n" >> provenance.yml
    printf -- "  tool_name: bcftools\\n  tool_version: \$(bcftools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> provenance.yml

    printf -- "- process_name: samtools\\n" >> provenance.yml
    printf -- "  tool_name: samtools\\n  tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> provenance.yml

	printf -- "- process_name: blast\\n" >> provenance.yml
	printf -- "  tool_name: blastn\\n  tool_version: \$(blastn -version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> provenance.yml

	printf -- "- process_name: spades\\n" >> provenance.yml
    printf -- "  tool_name: spades\\n  tool_version: \$(rnaviralspades.py --version 2>&1 | cut -d' ' -f4)\\n" >> provenance.yml

	"""
}

process pipeline_provenance {

	tag { pipeline_name + " / " + pipeline_version }

  	executor 'local'

	input:
	tuple val(pipeline_name), val(pipeline_version), val(analysis_start)

	output:
	file("pipeline_provenance.yml")

	script:
	"""
	printf -- "- pipeline_name: ${pipeline_name}\\n  pipeline_version: ${pipeline_version}\\n- timestamp_analysis_start: ${analysis_start}\\n" > pipeline_provenance.yml
	"""
}
