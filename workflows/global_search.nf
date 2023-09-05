#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { prep_database; make_blast_database; extract_genes_blast; run_self_blast; run_blastn; run_blastx; combine_references} from '../modules/blast.nf'
//include { p_make_blast_database; p_run_self_blast; p_run_blastn; p_filter_alignments; p_get_best_references; p_run_blastx } from './modules/p_blast.nf'
include { assembly } from '../modules/assembly.nf'
include { create_fasta_index; map_reads; sort_filter_index_sam; merge_fasta_bam} from '../modules/mapping.nf'
include { run_freebayes; run_mpileup ; get_common_snps } from '../modules/variant_calling.nf'
include { get_coverage; plot_coverage; make_pileup} from "../modules/coverage.nf"


workflow global_reference_search {
	take:
		ch_contigs
		ch_fastqs

	main:
		// DATABASE SETUP 
		ch_global_db_full = Channel.fromPath(params.reference_db_full)
		make_blast_database(ch_global_db_full).first().set{ch_global_db_full_blast}

		// SEARCH REFERENCE DB 
		run_blastn(ch_contigs.combine(ch_global_db_full), ch_global_db_full_blast, Channel.fromPath('NO_FILE').first())

		// SINGLE REFERENCE MAPPING 
		map_reads(ch_fastqs.join(run_blastn.out.ref))
		sort_filter_index_sam(map_reads.out)
		get_coverage(sort_filter_index_sam.out)

		// GENE DATABASE SETUP 
		extract_genes_blast(ch_global_db_full)
		

	emit:
		main = get_coverage.out.metric
			.join(run_blastn.out.ref)
			.join(sort_filter_index_sam.out)
			.join(get_coverage.out.main)
			.combine(ch_global_db_full.combine(extract_genes_blast.out.pos_vp1))
			.combine(ch_global_db_full.combine(extract_genes_blast.out.pos_rdrp))
			.map{ it -> [it[0], it.subList(1,it.size())]}
}
