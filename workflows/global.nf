#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { run_qualimap; run_custom_qc } from '../modules/qc.nf'
include { build_composite_reference; dehost_fastq } from '../modules/prep.nf'
include { prep_database; make_blast_database; extract_genes_database; run_self_blast; run_blastn; run_blastx} from '../modules/blast.nf'
include { combine_references} from '../modules/blast.nf'
include { assembly } from '../modules/assembly.nf'
include { map_reads; sort_filter_index_sam; merge_fasta_bam} from '../modules/mapping.nf'
include { get_coverage; plot_coverage; make_pileup} from "../modules/coverage.nf"
include { run_freebayes; run_mpileup ; get_common_snps } from '../modules/variant_calling.nf'
include { mask_low_coverage; make_consensus } from '../modules/consensus.nf'


workflow GLOBAL_ANALYSIS {
	take:
		ch_fastqs
		ch_contigs
		
	main:
		// DATABASE SETUP 
		ch_global_db_full = Channel.fromPath(params.global_database)
		prep_database(ch_global_db_full)
		make_blast_database(prep_database.out).first().set{ch_global_db_full_blast}

		// SEARCH REFERENCE DB 
		if (params.assemble) {
			run_blastn(ch_contigs.combine(ch_contigs), ch_global_db_full_blast, Channel.fromPath('NO_FILE').first())
		}else{
			run_blastn(ch_contigs.combine(prep_database.out), ch_global_db_full_blast, Channel.fromPath('NO_FILE').first())
		}

		blast_collect = run_blastn.out.filter.collectFile(name: "${params.outpath}/global_types_collect.tsv", keepHeader: true, skip: 1)

		// SINGLE REFERENCE MAPPING 
		map_reads(ch_fastqs.join(run_blastn.out.ref))
		sort_filter_index_sam(map_reads.out)

		ch_reference = run_blastn.out.ref
		ch_bamfile = sort_filter_index_sam.out

		// CALCULATE COVERAGE
		get_coverage(ch_bamfile)
		// PLOT COVERAGE
		plot_coverage(get_coverage.out.main)

		// COMPUTE PILEUP AND DETECT AMBIGUOUS POSITIONS
		make_pileup(ch_reference.join(ch_bamfile))

		// CHECK MAPPING QUALITY
		run_qualimap(ch_bamfile)

		// VARIANT CALLING
		run_freebayes(ch_bamfile.join(ch_reference))
		run_mpileup(ch_bamfile.join(ch_reference))
		get_common_snps(run_freebayes.out.join(run_mpileup.out))
		
		// CONSENSUS GENERATION 
		mask_low_coverage(ch_bamfile)
		make_consensus(get_common_snps.out.join(ch_reference).join(mask_low_coverage.out.bed))
	
		// CUSTOM QC 
		run_custom_qc(ch_bamfile.join(ch_reference).join(make_consensus.out.consensus))
		ch_qc_all = run_custom_qc.out.csv.collectFile(name: "${params.outpath}/qc/global_qc_all.csv", keepHeader: true, skip: 1)
		

	emit:
		// main = get_coverage.out.metric
		// 	.join(run_blastn.out.ref)
		// 	.join(sort_filter_index_sam.out)
		// 	.join(get_coverage.out.main)
		// 	.combine(ch_global_db_full.combine(extract_genes_database.out.pos_vp1))
		// 	.combine(ch_global_db_full.combine(extract_genes_database.out.pos_rdrp))
		// 	.map{ it -> [it[0], it.subList(1,it.size())]}
		blast_collect = blast_collect

}