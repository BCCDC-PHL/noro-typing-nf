#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { qualimap; custom_qc } from '../modules/qc.nf'
include { prep_database; make_blast_database; extract_genes_database; run_self_blast; blastn; combine_references} from '../modules/blast.nf'
include { assembly } from '../modules/assembly.nf'
include { map_reads; sort_filter_sam; merge_fasta_bam} from '../modules/mapping.nf'
include { get_coverage; plot_coverage; make_pileup} from "../modules/coverage.nf"
include { freebayes; mpileup ; find_common_snps } from '../modules/variant_calling.nf'
include { mask_low_coverage; make_consensus } from '../modules/consensus.nf'


workflow GLOBAL_ANALYSIS {
	take:
		ch_fastqs
		ch_contigs
		
	main:
		// DATABASE SETUP 
		ch_blastn_global_db = Channel.fromPath(params.blastn_global_database)
		prep_database(ch_blastn_global_db)
		make_blast_database(prep_database.out).first().set{ch_blastn_global_database}

		// SEARCH REFERENCE DB 
		if (params.assemble) {
			blastn(ch_contigs.combine(ch_contigs), ch_blastn_global_database, Channel.fromPath('NO_FILE').first())
		}else{
			blastn(ch_contigs.combine(prep_database.out), ch_blastn_global_database, Channel.fromPath('NO_FILE').first())
		}

		blast_collect = blastn.out.filter.collectFile(name: "${params.outdir}/global_types_collect.tsv", keepHeader: true, skip: 1)

		// SINGLE REFERENCE MAPPING 
		map_reads(ch_fastqs.join(blastn.out.ref))
		sort_filter_sam(map_reads.out.sam)

		ch_reference = blastn.out.ref
		ch_bamfile = sort_filter_sam.out.bam

		// CALCULATE COVERAGE
		get_coverage(ch_bamfile)
		// PLOT COVERAGE
		plot_coverage(get_coverage.out.main)

		// COMPUTE PILEUP AND DETECT AMBIGUOUS POSITIONS
		make_pileup(ch_reference.join(ch_bamfile))

		// CHECK MAPPING QUALITY
		qualimap(ch_bamfile)

		// VARIANT CALLING
		freebayes(ch_bamfile.join(ch_reference))
		mpileup(ch_bamfile.join(ch_reference))
		find_common_snps(freebayes.out.vcf.join(mpileup.out.vcf))
		
		// CONSENSUS GENERATION 
		mask_low_coverage(ch_bamfile)
		make_consensus(find_common_snps.out.vcf.join(ch_reference).join(mask_low_coverage.out.bed))
	
		// CUSTOM QC 
		custom_qc(ch_bamfile.join(ch_reference).join(make_consensus.out.consensus))
		ch_qc_all = custom_qc.out.csv.collectFile(name: "${params.outdir}/qc/global_qc_all.csv", keepHeader: true, skip: 1)
		
		ch_provenance = blastn.out.provenance
		ch_provenance = ch_provenance.mix(map_reads.out.provenance)
		ch_provenance = ch_provenance.mix(sort_filter_sam.out.provenance)
		//ch_provenance = ch_provenance.mix(get_coverage.out.provenance)
		ch_provenance = ch_provenance.mix(qualimap.out.provenance)
		ch_provenance = ch_provenance.mix(mpileup.out.provenance)
		ch_provenance = ch_provenance.mix(freebayes.out.provenance)
		ch_provenance = ch_provenance.mix(find_common_snps.out.provenance)

	emit:
		blast_collect = blast_collect
		provenance = ch_provenance

}