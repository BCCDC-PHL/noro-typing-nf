#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

// include { pipeline_provenance } from './modules/provenance.nf'
// include { collect_provenance } from './modules/provenance.nf'
include { fastqc; quast; make_typing_report } from './modules/qc.nf'
include { hash_files ; make_union_database; cutadapt; fastp; kraken2; filter_reads_kraken2 } from './modules/prep.nf' 
include { build_composite_reference; dehost_fastq } from './modules/prep.nf'
include { combine_references} from './modules/blast.nf'
include { assembly } from './modules/assembly.nf'
include { multiqc } from './modules/multiqc.nf'
include { collect_provenance } from './modules/provenance.nf'

include { COMPOSITE_ANALYSIS } from './workflows/composite.nf'
include { GLOBAL_ANALYSIS } from './workflows/global.nf'

println "HELLO. STARTING NOROVIRUS TYPING PIPELINE."
println "${LocalDateTime.now()}"



// prints to the screen and to the log
log.info """Norovirus Typing Pipeline
===================================
projectDir        : ${projectDir}
launchDir         : ${launchDir}
primers           : ${params.primers}
gtype_blast_db          : ${params.blastn_gtype_database}
ptype_blast_db          : ${params.blastn_ptype_database}
global_blast_db    : ${params.blastn_global_database}
fastqInputDir     : ${params.fastq_input}
outdir            : ${params.outdir}
pipeline_cache    : ${params.pipeline_cache}
run_name          : ${params.run_name}
git repo          : $workflow.repository
git version       : $workflow.revision [$workflow.commitId]
user              : $workflow.userName
""".stripIndent()

// database          : ${params.db}
// Git repository    : $workflow.repository
// git commit id     : $workflow.commitId
// branch            : $workflow.revision
// pipeline run      : ${params.pipeline_short_name}
// pipeline version  : ${params.pipeline_minor_version}



workflow {
	ch_start_time = Channel.of(LocalDateTime.now())
	ch_pipeline_name = Channel.of(workflow.manifest.name)
	ch_pipeline_version = Channel.of(workflow.manifest.version)

	//ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

	ch_human_ref = Channel.from(params.human_reference_fasta)
	ch_blastn_gtype_database = Channel.from(params.blastn_gtype_database)
	ch_blastn_ptype_database = Channel.from(params.blastn_ptype_database)
	ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }

	main:
		hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))

		if ("bsr" in params.blast_metrics_global.split(",")){
			error "ERROR: Global reference search cannot use BSR (blast score ratio)."
		}

		ch_sample_list = ch_fastq_input.map{it -> it[0]}.collectFile(name: "${params.outdir}/sample-list.txt", newLine: true)
		
		// READ TRIMMING AND FILTERING 
		fastp(ch_fastq_input)
		cutadapt(fastp.out.trimmed_reads)
		
		// FASTQ QUALITY CONTROL 
		fastqc(fastp.out.trimmed_reads)
		
		// RUN KRAKEN2 WITH OPTIONAL READ FILTERING
		kraken2(fastp.out.trimmed_reads)

		if (params.k2_read_filter) {
			ch_filtered_reads = filter_reads_kraken2(fastp.out.trimmed_reads.join(kraken2.main)).out.fastq
		}else{
			ch_filtered_reads = fastp.out.trimmed_reads
		}

		
		// DEHOSTING
		make_union_database(ch_blastn_gtype_database, ch_blastn_ptype_database).first()
		build_composite_reference(ch_human_ref.combine(make_union_database.out.fasta))
		dehost_fastq(
			ch_filtered_reads,
			make_union_database.out.headers.first(), 
			build_composite_reference.out.first()
		) 

		// ASSEMBLY & QC
		assembly(dehost_fastq.out.fastq)
		quast(assembly.out.contigs.map{ it -> it[1]}.collect())

		// Run analysis using composite reference (competitive mapping to two references simultaneously)
		COMPOSITE_ANALYSIS(dehost_fastq.out.fastq, assembly.out.contigs)

		// If global database is passed, use full length global reference database to find reference
		if (params.blastn_global_database){
			// perform a BLAST search on the global reference database 
			GLOBAL_ANALYSIS(dehost_fastq.out.fastq, assembly.out.contigs)
			ch_global_analysis = GLOBAL_ANALYSIS.out.blast_collect

		} else {
			ch_global_analysis = Channel.from('NO_FILE').first()
		}


		// Make a master output report 
		make_typing_report(
			ch_sample_list,
			COMPOSITE_ANALYSIS.out.qc_all,
			COMPOSITE_ANALYSIS.out.gblast_collect,
			COMPOSITE_ANALYSIS.out.pblast_collect,
			ch_global_analysis
		)

		// Collect all the relevant files for MULTIQC
		ch_multiqc_inputs = Channel.empty()
		ch_multiqc_inputs = ch_multiqc_inputs.mix(fastp.out.json.map{it -> it[1]})
		ch_multiqc_inputs = ch_multiqc_inputs.mix(fastqc.out.zip.map{ it -> [it[1], it[2]]}.collect())
		ch_multiqc_inputs = ch_multiqc_inputs.mix(cutadapt.out.log)
		ch_multiqc_inputs = ch_multiqc_inputs.mix(kraken2.out.report)
		ch_multiqc_inputs = ch_multiqc_inputs.mix(quast.out.tsv)
		//ch_multiqc_inputs = ch_multiqc_inputs.mix(qualimap.out.main)
		multiqc(ch_multiqc_inputs.collect().ifEmpty([]) )

		ch_provenance = hash_files.out.provenance
		ch_provenance = ch_provenance.mix(cutadapt.out.provenance)
		ch_provenance = ch_provenance.mix(fastp.out.provenance)
		ch_provenance = ch_provenance.mix(fastqc.out.provenance)
		ch_provenance = ch_provenance.mix(kraken2.out.provenance)
		ch_provenance = ch_provenance.mix(dehost_fastq.out.provenance)
		ch_provenance = ch_provenance.mix(assembly.out.provenance)
		ch_provenance = ch_provenance.mix(COMPOSITE_ANALYSIS.out.provenance)
		ch_provenance = ch_provenance.mix(GLOBAL_ANALYSIS.out.provenance).unique().groupTuple()
		ch_provenance.view()
		collect_provenance(ch_provenance)
}