#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { build_composite_reference; dehost_fastq } from './modules/prep.nf'
include { prep_database; make_blast_database; extract_genes_blast; run_self_blast; run_blastn; run_blastx} from './modules/blast.nf'


workflow genotyping {

	take:
		ch_contigs
		ch_blastdb_fasta
	main:
		// GENOTYPE BLAST SEARCH
		prep_database(ch_blastdb_fasta)
		extract_genes_blast(prep_database.out)
		make_blast_database(extract_genes_blast.out.db).first().set{ch_blastdb} // first() converts the channel from a consumable queue channel into an infinite single-value channel
		// run_blastn(ch_contigs, ch_blastdb)

		// Needed for blast score ratio
		if (params.blast_typing_metric == 'bsr') {
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
		gene_db = extract_genes_blast.out.db
		db_pos = prep_database.out.combine(extract_genes_blast.out.pos)
}

workflow ptyping {
	take:
		ch_contigs
		ch_blastdb_fasta
	main:
		// GENOTYPE BLAST SEARCH
		prep_database(ch_blastdb_fasta)
		extract_genes_blast(prep_database.out)
		make_blast_database(extract_genes_blast.out.db).first().set{ch_blastdb} // first() converts the channel from a consumable queue channel into an infinite single-value channel

		// Needed for blast score ratio
		if (params.blast_typing_metric == 'bsr') {
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
		gene_db = extract_genes_blast.out.db
		db_pos = prep_database.out.combine(extract_genes_blast.out.pos)
}