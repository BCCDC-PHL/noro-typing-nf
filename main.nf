#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

// include { pipeline_provenance } from './modules/provenance.nf'
// include { collect_provenance } from './modules/provenance.nf'
include { fastQC; fastq_check; run_quast; run_qualimap} from './modules/qc.nf'
include { merge_databases; cutadapt; fastp; fastp_json_to_csv; run_kraken; kraken_filter } from './modules/prep.nf' 
include { build_composite_reference; index_composite_reference ; get_reference_headers; dehost_fastq } from './modules/prep.nf'
include { make_blast_database; run_self_blast; run_blastn; filter_alignments; run_blastx; select_best_reference } from './modules/blast.nf'
//include { p_make_blast_database; p_run_self_blast; p_run_blastn; p_filter_alignments; p_get_best_references; p_run_blastx } from './modules/p_blast.nf'
include { assembly } from './modules/assembly.nf'
include { create_bwa_index; create_fasta_index; map_reads; sort_filter_index_sam; index_bam } from './modules/mapping.nf'
include { run_freebayes; run_mpileup ; get_common_snps } from './modules/variant_calling.nf'
include { get_coverage; plot_coverage} from "./modules/coverage.nf"
include { mask_low_coverage; make_consensus } from './modules/consensus.nf'
include { make_multifasta; make_msa; make_tree } from './modules/phylo.nf'
include { multiqc } from './modules/multiqc.nf'

println "HELLO. STARTING NOROVIRUS METAGENOMICS PIPELINE."
println "${LocalDateTime.now()}"

// prints to the screen and to the log
log.info """Norovirus Metagenomics Pipeline
===================================
projectDir        : ${projectDir}
launchDir         : ${launchDir}
primers           : ${params.primers}
blast_db          : ${params.blastdb_gtype_fasta}
fastqInputDir     : ${params.fastq_input}
outdir            : ${params.outdir}
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

blast_db_dir = file("${projectDir}/cache/blast_db").mkdirs()
composite_dir = file("${projectDir}/cache/composite").mkdirs()
println blast_db_dir ? "BLAST DB directory created successfully" : "Cannot create directory: $blast_db_dir"
println composite_dir ? "Composite directory created successfully" : "Cannot create directory: $composite_dir"

// println "${params.blastdb_gtype_fasta}"
// println "${params.blastdb_gtype_fasta}" =~ /blast/ ? "SUCCESS" : "NOT FOUND" 

workflow genotyping {

	take:
		ch_contigs
		ch_blastdb_fasta
	main:
		// GENOTYPE BLAST SEARCH
		def workflow_type = 'gtype'
		make_blast_database(ch_blastdb_fasta).first().set{ch_blastdb} // first() converts the channel from a consumable queue channel into an infinite single-value channel
		run_blastn(ch_contigs, ch_blastdb)

		// Needed for blast score ratio
		run_self_blast(ch_blastdb)

		if (params.assemble){
			filter_alignments(run_blastn.out.join(ch_contigs).combine(run_self_blast.out))
		} else {
			filter_alignments(run_blastn.out.combine(ch_blastdb_fasta).combine(run_self_blast.out))
		}

		//get_best_references(ch_gtype_refs)
	emit:
		main = filter_alignments.out.main
		filter = filter_alignments.out.filter
}

workflow ptyping {
	take:
		ch_contigs
		ch_blastdb_fasta
	main:
		// GENOTYPE BLAST SEARCH
		params.workflow = 'ptype'
		make_blast_database(ch_blastdb_fasta).first().set{ch_blastdb} // first() converts the channel from a consumable queue channel into an infinite single-value channel
		run_blastn(ch_contigs, ch_blastdb)

		// Needed for blast score ratio
		run_self_blast(ch_blastdb)

		if (params.assemble){
			filter_alignments(run_blastn.out.join(ch_contigs).combine(run_self_blast.out))
		} else {
			filter_alignments(run_blastn.out.combine(ch_blastdb_fasta).combine(run_self_blast.out))
		}
		
		//get_best_references(filter_alignments.out)
	emit:
		main = filter_alignments.out.main
		filter = filter_alignments.out.filter
}

workflow {
	ch_start_time = Channel.of(LocalDateTime.now())
	ch_pipeline_name = Channel.of(workflow.manifest.name)
	ch_pipeline_version = Channel.of(workflow.manifest.version)

	//ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

	//ch_ref_names = Channel.fromList(params.virus_ref_names).collect()
	ch_composite_paths = Channel.fromList([params.human_ref, params.virus_ref]).collect()
	ch_human_ref = Channel.from(params.human_ref)
	ch_centrifuge_db = Channel.from(params.centrifuge_db)
	ch_blastdb_gtype_fasta = Channel.from(params.blastdb_gtype_fasta)
	ch_blastdb_ptype_fasta = Channel.from(params.blastdb_ptype_fasta)
	ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }

	main:
		//hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))
		
		// UNION DATABASE
		merge_databases(ch_blastdb_gtype_fasta, ch_blastdb_ptype_fasta).set{ch_union_db}

		// QUALITY CONTROL 
		cutadapt(ch_fastq_input)
		fastp(cutadapt.out.trimmed_reads)
		// fastp(ch_fastq_input)
		fastp_json_to_csv(fastp.out.json)
		
		fastQC(fastp.out.trimmed_reads)
		fastq_check(fastp.out.trimmed_reads)
		fastq_check.out.formatted.collectFile(name: "${params.outdir}/fastq_qual/fastq_stats.tsv", keepHeader: true, skip: 1)
		
		// KRAKEN FILTERING
		run_kraken(fastp.out.trimmed_reads)
		kraken_filter(fastp.out.trimmed_reads.join(run_kraken.out))
		
		// DEHOSTING
		build_composite_reference(ch_human_ref.combine(ch_union_db))
		index_composite_reference(build_composite_reference.out)
		get_reference_headers(ch_union_db)
		dehost_fastq(
			kraken_filter.out,
			get_reference_headers.out.first(), 
			index_composite_reference.out.first()
		) 

		// ASSEMBLY
		assembly(dehost_fastq.out.fastq)
		run_quast(assembly.out.map{ it -> it[1]}.collect())

		// GENOTYPING / PTYPING 
		genotyping(assembly.out, ch_blastdb_gtype_fasta)
		ptyping(assembly.out, ch_blastdb_ptype_fasta)

		// Collect typing results 
		genotyping.out.filter.collectFile(name: "${params.outdir}/blastn/final/gtypes.tsv", keepHeader: true, skip: 1)
		ptyping.out.filter.collectFile(name: "${params.outdir}/blastn/final/ptypes.tsv", keepHeader: true, skip: 1)

		// Pick best reference out of G & P type candidates
		select_best_reference(genotyping.out.main.join(ptyping.out.main))

		select_best_reference.out.blast.collect().collectFile(name: "${params.outdir}/blastn/final/final_types.tsv", keepHeader: true, skip: 1).set{ch_types}
		
		create_bwa_index(select_best_reference.out.ref)

		// READ MAPPING, SORTING, FILTERING
		map_reads(fastp.out.trimmed_reads.join(create_bwa_index.out))
		sort_filter_index_sam(map_reads.out)

		run_qualimap(sort_filter_index_sam.out)

		// // PLOT COVERAGE
		// get_coverage(sort_filter_index_sam.out)
		// plot_coverage(get_coverage.out.coverage_file)

		// // VARIANT CALLING
		// create_fasta_index(select_best_reference.out.ref)
		// run_freebayes(sort_filter_index_sam.out.join(create_fasta_index.out))
		// run_mpileup(sort_filter_index_sam.out.join(create_fasta_index.out))
		// get_common_snps(run_freebayes.out.join(run_mpileup.out))
		
		// // CONSENSUS GENERATION 
		// mask_low_coverage(sort_filter_index_sam.out)
		// make_consensus(get_common_snps.out.join(select_best_reference.out.ref).join(mask_low_coverage.out))

		// make_multifasta(make_consensus.out.collect())
		// make_msa(make_multifasta.out)
		// make_tree(make_msa.out)


		// QualiMap(FluViewer.out.alignment)
		// parseQMresults(QualiMap.out.genome_results)

		// Collect al the relevant filesfor MULTIQC
		ch_multiqc_inputs = Channel.empty()
		
		ch_multiqc_inputs.mix(fastQC.out.zip.map{ it -> [it[1], it[2]]}.collect()).set{ch_multiqc_inputs}
		ch_multiqc_inputs.mix(cutadapt.out.log).set{ch_multiqc_inputs}
		ch_multiqc_inputs.mix(run_quast.out.tsv).set{ch_multiqc_inputs}
		ch_multiqc_inputs.mix(fastp.out.json.map{it -> it[1]}).set{ch_multiqc_inputs}
		ch_multiqc_inputs.mix(run_qualimap.out.main).set{ch_multiqc_inputs}
		multiqc(ch_multiqc_inputs.collect().ifEmpty([]) )
		
		// segcov(FluViewer.out.alignment)

		// ch_provenance = FluViewer.out.provenance
		// ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], [it[1]] << it[2]] }
		// ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(cutadapt.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(ch_fastq_input.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
		// collect_provenance(ch_provenance)

}