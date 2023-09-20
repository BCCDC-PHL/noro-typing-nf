#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { build_composite_reference; dehost_fastq } from '../modules/prep.nf'
include { prep_database; make_blast_database; extract_genes_blast; run_self_blast; run_blastn; run_blastx} from '../modules/blast.nf'
include { combine_references} from '../modules/blast.nf'
include { assembly } from '../modules/assembly.nf'
include { create_fasta_index; map_reads; sort_filter_index_sam; merge_fasta_bam} from '../modules/mapping.nf'
include { get_coverage; plot_coverage; make_pileup} from "../modules/coverage.nf"

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

		// Collect typing results 
		blast_collect = run_blastn.out.filter.collectFile(name: "${params.outdir}/blastn/gtype/filtered/gtypes.tsv", keepHeader: true, skip: 1)

		// ch_gtype_blast_results = Channel.from("${params.outdir}/blastn/gtype/filtered/gtypes.tsv").first()


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
		blast_collect = run_blastn.out.filter.collectFile(name: "${params.outdir}/blastn/ptype/filtered/ptypes.tsv", keepHeader: true, skip: 1)
}


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
		blast_collect = run_blastn.out.filter.collectFile(name: "${params.outdir}/blastn/global/filtered/global_types.tsv", keepHeader: true, skip: 1)

}
