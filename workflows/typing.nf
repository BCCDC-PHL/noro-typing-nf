#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { prep_database; make_blast_database; extract_genes_database; run_self_blast; blastn} from '../modules/blast.nf'

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
			blastn(ch_contigs.join(ch_contigs), ch_blastdb, ch_selfblast)
		} else {
			blastn(ch_contigs.combine(prep_database.out), ch_blastdb, ch_selfblast)
		}

	emit:
		main = blastn.out.main
		filter = blastn.out.filter
		provenance = blastn.out.provenance

}

