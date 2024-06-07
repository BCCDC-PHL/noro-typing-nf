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

include { TYPING as G_TYPING } from "./typing.nf" 
include { TYPING as P_TYPING } from "./typing.nf" 


workflow COMPOSITE_ANALYSIS {

	take:
		ch_fastq
		ch_contigs

	main:
		ch_blastdb_gtype_fasta = Channel.from(params.gtype_database)
		ch_blastdb_ptype_fasta = Channel.from(params.ptype_database)

		// GENOTYPING / PTYPING 
		G_TYPING(ch_contigs, ch_blastdb_gtype_fasta)
		P_TYPING(ch_contigs, ch_blastdb_ptype_fasta)

		gblast_collect = G_TYPING.out.filter.collectFile(name: "${params.outpath}/composite_gtypes_collect.tsv", keepHeader: true, skip: 1)
		pblast_collect = P_TYPING.out.filter.collectFile(name: "${params.outpath}/composite_ptypes_collect.tsv", keepHeader: true, skip: 1)


		// Create combined reference
		// COMBINE REFERENCES
		combine_references(G_TYPING.out.main.join(P_TYPING.out.main))

		// COMPETITIVE MAPPING 
		map_reads(ch_fastq.join(combine_references.out))
		sort_filter_index_sam(map_reads.out)

		// GENERATE NEW MERGED REFERENCE
		merge_fasta_bam(combine_references.out.join(sort_filter_index_sam.out))

		ch_best_reference = merge_fasta_bam.out.ref
		ch_bamfile = merge_fasta_bam.out.bam

		// CALCULATE COVERAGE
		get_coverage(ch_bamfile)
		// PLOT COVERAGE
		plot_coverage(get_coverage.out.main)

		// COMPUTE PILEUP AND DETECT AMBIGUOUS POSITIONS
		make_pileup(ch_best_reference.join(ch_bamfile))

		// CHECK MAPPING QUALITY
		run_qualimap(ch_bamfile)

		// VARIANT CALLING
		run_freebayes(ch_bamfile.join(ch_best_reference))
		run_mpileup(ch_bamfile.join(ch_best_reference))
		get_common_snps(run_freebayes.out.join(run_mpileup.out))
		
		// CONSENSUS GENERATION 
		mask_low_coverage(ch_bamfile)
		make_consensus(get_common_snps.out.join(ch_best_reference).join(mask_low_coverage.out.bed))
	
		// CUSTOM QC 
		run_custom_qc(ch_bamfile.join(ch_best_reference).join(make_consensus.out.consensus))
		ch_qc_all = run_custom_qc.out.csv.collectFile(name: "${params.outpath}/qc/composite_qc_all.csv", keepHeader: true, skip: 1)

		// ch_provenance = FluViewer.out.provenance
		// ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], [it[1]] << it[2]] }
		// ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(cutadapt.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(ch_fastq_input.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
		// collect_provenance(ch_provenance)

	emit:
		qc_all = ch_qc_all
		gblast_collect = gblast_collect
		pblast_collect = pblast_collect

}
