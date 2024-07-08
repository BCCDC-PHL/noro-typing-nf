#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { qualimap; custom_qc } from '../modules/qc.nf'
include { prep_database; make_blast_database; combine_references; extract_genes_database; run_self_blast; blastn} from '../modules/blast.nf'
include { assembly } from '../modules/assembly.nf'
include { map_reads; sort_filter_sam; merge_fasta_bam} from '../modules/mapping.nf'
include { get_coverage; plot_coverage; make_pileup} from "../modules/coverage.nf"
include { freebayes; mpileup ; find_common_snps } from '../modules/variant_calling.nf'
include { mask_low_coverage; make_consensus } from '../modules/consensus.nf'

include { TYPING as G_TYPING } from "./typing.nf" 
include { TYPING as P_TYPING } from "./typing.nf" 


workflow COMPOSITE_ANALYSIS {

	take:
		ch_fastq
		ch_contigs

	main:
		ch_blastn_gtype_database = Channel.from(params.blastn_gtype_database)
		ch_blastn_ptype_database = Channel.from(params.blastn_ptype_database)

		// GENOTYPING / PTYPING 
		G_TYPING(ch_contigs, ch_blastn_gtype_database)
		P_TYPING(ch_contigs, ch_blastn_ptype_database)

		gblast_collect = G_TYPING.out.filter.collectFile(name: "${params.outdir}/composite_gtypes_collect.tsv", keepHeader: true, skip: 1)
		pblast_collect = P_TYPING.out.filter.collectFile(name: "${params.outdir}/composite_ptypes_collect.tsv", keepHeader: true, skip: 1)


		// Create combined reference
		// COMBINE REFERENCES
		combine_references(G_TYPING.out.main.join(P_TYPING.out.main))

		// COMPETITIVE MAPPING 
		map_reads(ch_fastq.join(combine_references.out))
		sort_filter_sam(map_reads.out.sam)

		// GENERATE NEW MERGED REFERENCE
		merge_fasta_bam(combine_references.out.join(sort_filter_sam.out.bam))

		ch_best_reference = merge_fasta_bam.out.ref
		ch_bamfile = merge_fasta_bam.out.bam

		// CALCULATE COVERAGE
		get_coverage(ch_bamfile)
		// PLOT COVERAGE
		plot_coverage(get_coverage.out.main)

		// COMPUTE PILEUP AND DETECT AMBIGUOUS POSITIONS
		make_pileup(ch_best_reference.join(ch_bamfile))

		// CHECK MAPPING QUALITY
		qualimap(ch_bamfile)

		// VARIANT CALLING
		mpileup(ch_bamfile.join(ch_best_reference))
		freebayes(ch_bamfile.join(ch_best_reference))
		find_common_snps(freebayes.out.vcf.join(mpileup.out.vcf))
		
		// CONSENSUS GENERATION 
		mask_low_coverage(ch_bamfile)
		make_consensus(find_common_snps.out.vcf.join(ch_best_reference).join(mask_low_coverage.out.bed))
	
		// CUSTOM QC 
		custom_qc(ch_bamfile.join(ch_best_reference).join(make_consensus.out.consensus))
		ch_qc_all = custom_qc.out.csv.collectFile(name: "${params.outdir}/qc/composite_qc_all.csv", keepHeader: true, skip: 1)

		ch_provenance = G_TYPING.out.provenance
		ch_provenance = ch_provenance.mix(P_TYPING.out.provenance)
		ch_provenance = ch_provenance.mix(map_reads.out.provenance)
		ch_provenance = ch_provenance.mix(sort_filter_sam.out.provenance)
		//ch_provenance = ch_provenance.mix(get_coverage.out.provenance)
		ch_provenance = ch_provenance.mix(qualimap.out.provenance)
		ch_provenance = ch_provenance.mix(mpileup.out.provenance)
		ch_provenance = ch_provenance.mix(freebayes.out.provenance)
		ch_provenance = ch_provenance.mix(find_common_snps.out.provenance)


	emit:
		qc_all = ch_qc_all
		gblast_collect = gblast_collect
		pblast_collect = pblast_collect
		provenance = ch_provenance

}
