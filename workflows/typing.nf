#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { build_composite_reference; dehost_fastq } from '../modules/prep.nf'
include { prep_database; make_blast_database; extract_genes_database; run_self_blast; run_blastn; run_blastx} from '../modules/blast.nf'
include { combine_references} from '../modules/blast.nf'
include { assembly } from '../modules/assembly.nf'
include { map_reads; sort_filter_index_sam; merge_fasta_bam} from '../modules/mapping.nf'
include { get_coverage; plot_coverage; make_pileup} from "../modules/coverage.nf"

workflow TYPING {

	take:
		ch_contigs
		ch_blastdb_fasta
	main:
		// GENOTYPE BLAST SEARCH
		prep_database(ch_blastdb_fasta)
		extract_genes_database(prep_database.out)
		make_blast_database(extract_genes_database.out.db).first().set{ch_blastdb} // first() converts the channel from a consumable queue channel into an infinite single-value channel

		// Needed for blast score ratio
		if (params.blast_metrics_composite == 'bsr') {
			ch_selfblast = run_self_blast(ch_blastdb)
		} else {
			ch_selfblast = Channel.fromPath("NO_FILE").first()
		}

		if (params.assemble){
			run_blastn(ch_contigs.join(ch_contigs), ch_blastdb, ch_selfblast)
		} else {
			run_blastn(ch_contigs.combine(prep_database.out), ch_blastdb, ch_selfblast)
		}

	emit:
		main = run_blastn.out.main
		filter = run_blastn.out.filter

}

