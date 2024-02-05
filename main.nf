#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

// include { pipeline_provenance } from './modules/provenance.nf'
// include { collect_provenance } from './modules/provenance.nf'
include { fastQC; fastq_check; run_quast} from './modules/qc.nf'
include { merge_databases; cutadapt; fastp; fastp_json_to_csv; run_kraken; kraken_filter; run_centrifuge } from './modules/prep.nf' 
include { build_composite_reference; index_composite_reference ; get_reference_headers; dehost_fastq } from './modules/prep.nf'
include { make_blast_database; run_self_blast; run_blastn; run_blastx; filter_alignments; get_best_references; } from './modules/blast.nf'
//include { p_make_blast_database; p_run_self_blast; p_run_blastn; p_filter_alignments; p_get_best_references; p_run_blastx } from './modules/p_blast.nf'
include { run_spades; run_shovill; run_concoct; run_metabat } from './modules/assembly.nf'
include { create_bwa_index; create_fasta_index; map_reads; sort_filter_index_sam; index_bam } from './modules/mapping.nf'
include { run_freebayes; run_mpileup ; get_common_snps } from './modules/variant_calling.nf'
include { mask_low_coverage; make_consensus } from './modules/consensus.nf'
include { get_coverage; plot_coverage } from './modules/coverage.nf'
include { make_multifasta; make_msa; make_tree } from './modules/phylo.nf'
include { multiqc } from './modules/multiqc.nf'

println "HELLO. STARTING NOROVIRUS METAGENOMICS PIPELINE."
println "${LocalDateTime.now()}"

// prints to the screen and to the log
log.info """Norovirus Metagenomics Pipeline
===================================
projectDir        : ${projectDir}
launchDir         : ${launchDir}
fastqInputDir     : ${params.fastq_input}
outdir            : ${params.outdir}
run_name          : ${params.run_name}
git repo          : $workflow.repository
git version       : $workflow.revision [$workflow.commitId]
user              : $workflow.userName
""".stripIndent()


workflow {
	ch_start_time = Channel.of(LocalDateTime.now())
	ch_pipeline_name = Channel.of(workflow.manifest.name)
	ch_pipeline_version = Channel.of(workflow.manifest.version)

	//ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

	ch_adapters = Channel.fromPath(params.adapters_path)
	//ch_ref_names = Channel.fromList(params.virus_ref_names).collect()
	ch_centrifuge_db = Channel.from(params.centrifuge_db)
	ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }

	main:
		//hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))
		
		// UNION DATABASE
		//merge_databases(ch_blastdb_gtype_fasta, ch_blastdb_ptype_fasta).set{ch_union_db}
		//ch_fastq_input.combine(ch_adapters).subscribe{println "value: $it"}
		// QUALITY CONTROL 
		cutadapt(ch_fastq_input.combine(ch_adapters))
		fastp(cutadapt.out.trimmed_reads)
		// fastp(ch_fastq_input)
		fastp_json_to_csv(fastp.out.json)
		
		fastQC(fastp.out.trimmed_reads)
		fastq_check(fastp.out.trimmed_reads)
		fastq_check.out.formatted.collectFile(name: "${params.outdir}/fastq_qual/fastq_stats.tsv", keepHeader: true, skip: 1)
		
		run_shovill(fastp.out.trimmed_reads)
		create_bwa_index(run_shovill.out)
		map_reads(fastp.out.trimmed_reads.join(create_bwa_index.out))
		sort_filter_index_sam(map_reads.out)
		run_concoct(run_shovill.out.join(sort_filter_index_sam.out))
		run_metabat(run_shovill.out.join(sort_filter_index_sam.out))

		get_coverage(sort_filter_index_sam.out)
		plot_coverage(get_coverage.out.coverage_file)

		// // KRAKEN FILTERING
		run_kraken(fastp.out.trimmed_reads)
		run_centrifuge(fastp.out.trimmed_reads)
		// run_blastx(run_shovill.out)
		// kraken_filter(fastp.out.trimmed_reads.join(run_kraken.out))
		

		// genotyping( assembly.out, ch_blastdb_gtype_fasta)
		// ptyping(assembly.out, ch_blastdb_ptype_fasta)

		// create_bwa_index(genotyping.out.refs)

		

		// // READ MAPPING, SORTING, FILTERING
		// map_reads(fastp.out.trimmed_reads.join(create_bwa_index.out))
		// sort_filter_sam(map_reads.out)
		// index_bam(sort_filter_sam.out)

		// // VARIANT CALLING
		// create_fasta_index(genotyping.out.refs)
		// run_freebayes(sort_filter_sam.out.join(create_fasta_index.out))
		// run_mpileup(sort_filter_sam.out.join(create_fasta_index.out))
		// get_common_snps(run_freebayes.out.join(run_mpileup.out))
		
		// // CONSENSUS GENERATION 
		// mask_low_coverage(sort_filter_sam.out)
		// make_consensus(get_common_snps.out.join(genotyping.out.refs).join(mask_low_coverage.out))

		// make_multifasta(make_consensus.out.collect())
		// make_msa(make_multifasta.out)
		// make_tree(make_msa.out)


		// QualiMap(FluViewer.out.alignment)
		// parseQMresults(QualiMap.out.genome_results)

		// Collect al the relevant filesfor MULTIQC
		// ch_fastqc_collected = fastQC.out.zip.map{ it -> [it[1], it[2]]}.collect()
		// multiqc(fastp.out.json.mix( cutadapt.out.log, ch_fastqc_collected ).collect().ifEmpty([]) )

		// ch_provenance = FluViewer.out.provenance
		// ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], [it[1]] << it[2]] }
		// ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(cutadapt.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(ch_fastq_input.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
		// collect_provenance(ch_provenance)

}
