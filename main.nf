#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

// include { pipeline_provenance } from './modules/provenance.nf'
// include { collect_provenance } from './modules/provenance.nf'
include { fastp ; cutadapt; fastQC; fastq_check; } from './modules/qc.nf'
include { make_blast_database; run_blastn} from './modules/blast.nf'
include { assembly } from './modules/assembly.nf'
include { multiqc } from './modules/multiqc.nf'

println "HELLO. STARTING NOROVIRUS METAGENOMICS PIPELINE."

// prints to the screen and to the log
log.info """Norovirus Metagenomics Pipeline
		===================================
		projectDir        : ${projectDir}
		launchDir         : ${launchDir}
		primers           : ${params.primers}
		blast_db          : ${params.blast_db}
		fastqInputDir     : ${params.fastq_input}
		outdir            : ${params.outdir}
		run_name          : ${params.run_name}
		user              : $workflow.userName
		""".stripIndent()

// database          : ${params.db}
// Git repository    : $workflow.repository
// git commit id     : $workflow.commitId
// branch            : $workflow.revision
// pipeline run      : ${params.pipeline_short_name}
// pipeline version  : ${params.pipeline_minor_version}

blast_db_dir = file("${projectDir}/blast_db").mkdirs()
println blast_db_dir ? "BLAST DB directory created successfully" : "Cannot create directory: $blast_db_dir"

workflow {
	ch_start_time = Channel.of(LocalDateTime.now())
	ch_pipeline_name = Channel.of(workflow.manifest.name)
	ch_pipeline_version = Channel.of(workflow.manifest.version)

	//ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

	ch_primers = Channel.fromPath(params.adapters_path)
	ch_blast_db = Channel.fromPath(params.blast_db)
	ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }

	main:
		//hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))
		
		// QUALITY CONTROL 
		fastp( ch_fastq_input )
		// cutadapt(fastp.out.trimmed_reads.combine(ch_primers))
		fastQC(fastp.out.trimmed_reads)
		fastq_check(fastp.out.trimmed_reads)
		fastq_check.out.formatted.collectFile(name: "${params.outdir}/fastq_qual/fastq_stats.tsv", keepHeader: true, skip: 1)
		kraken(fastp.out.trimmed_reads)

		// ASSEMBLY
		assembly(fastp.out.trimmed_reads)

		// BLAST SEARCH
		make_blast_database(ch_blast_db)
		run_blastn(assembly.out.combine(make_blast_database.out))
		
		// FluViewer(cutadapt.out.primer_trimmed_reads.combine(ch_db))
		// FluViewer_assemble(cutadapt.out.primer_trimmed_reads.combine(ch_db))
		// QualiMap(FluViewer.out.alignment)
		// parseQMresults(QualiMap.out.genome_results)

		// Collect al the relevant filesfor MULTIQC
		// ch_fastqc_collected = fastQC.out.zip.map{ it -> [it[1], it[2]]}.collect()
		// multiqc(fastp.out.json.mix( cutadapt.out.log, ch_fastqc_collected ).collect().ifEmpty([]) )
		
		// segcov(FluViewer.out.alignment)

		// ch_provenance = FluViewer.out.provenance
		// ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], [it[1]] << it[2]] }
		// ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(cutadapt.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(ch_fastq_input.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
		// collect_provenance(ch_provenance)

}