#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

// include { pipeline_provenance } from './modules/provenance.nf'
// include { collect_provenance } from './modules/provenance.nf'
include { fastQC; fastq_check; run_quast; run_qualimap; run_custom_qc} from './modules/qc.nf'
include { make_union_database; cutadapt; fastp; fastp_json_to_csv; run_kraken; kraken_filter } from './modules/prep.nf' 
include { build_composite_reference; index_composite_reference ; dehost_fastq } from './modules/prep.nf'
include { prep_database; make_blast_database; extract_genes_blast; run_self_blast; run_blastn; run_blastx; combine_references ; find_reference} from './modules/blast.nf'
//include { p_make_blast_database; p_run_self_blast; p_run_blastn; p_filter_alignments; p_get_best_references; p_run_blastx } from './modules/p_blast.nf'
include { assembly } from './modules/assembly.nf'
include { create_bwa_index; create_fasta_index; map_reads; sort_filter_index_sam; merge_fasta_bam} from './modules/mapping.nf'
include { run_freebayes; run_mpileup ; get_common_snps } from './modules/variant_calling.nf'
include { get_coverage; plot_coverage; make_pileup} from "./modules/coverage.nf"
include { mask_low_coverage; make_consensus } from './modules/consensus.nf'
include { make_multifasta; make_msa; make_dates_file; make_tree; extract_sample_genes} from './modules/phylo.nf'
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
		prep_database(ch_blastdb_fasta)
		extract_genes_blast(prep_database.out)
		make_blast_database(extract_genes_blast.out).first().set{ch_blastdb} // first() converts the channel from a consumable queue channel into an infinite single-value channel
		// run_blastn(ch_contigs, ch_blastdb)

		// Needed for blast score ratio
		run_self_blast(ch_blastdb)

		if (params.assemble){
			run_blastn(ch_contigs.join(ch_contigs).combine(run_self_blast.out), ch_blastdb)
		} else {
			run_blastn(ch_contigs.combine(prep_database.out).combine(run_self_blast.out), ch_blastdb)
		}

	emit:
		main = run_blastn.out.main
		filter = run_blastn.out.filter
		gene_db = extract_genes_blast.out
}

workflow ptyping {
	take:
		ch_contigs
		ch_blastdb_fasta
	main:
		// GENOTYPE BLAST SEARCH
		prep_database(ch_blastdb_fasta)
		extract_genes_blast(prep_database.out)
		make_blast_database(extract_genes_blast.out).first().set{ch_blastdb} // first() converts the channel from a consumable queue channel into an infinite single-value channel

		// Needed for blast score ratio
		run_self_blast(ch_blastdb)

		if (params.assemble){
			run_blastn(ch_contigs.join(ch_contigs).combine(run_self_blast.out), ch_blastdb)
		} else {
			run_blastn(ch_contigs.combine(prep_database.out).combine(run_self_blast.out), ch_blastdb)
		}

	emit:
		main = run_blastn.out.main
		filter = run_blastn.out.filter
		gene_db = extract_genes_blast.out
}

workflow create_gtree {
	take:
		ch_consensus_fasta
		ch_gene_database_fasta

	main:
		extract_sample_genes(ch_consensus_fasta.combine(ch_gene_database_fasta))
		make_multifasta(ch_gene_database_fasta.mix(extract_sample_genes.out).collect())
		make_dates_file(make_multifasta.out)
		make_msa(make_multifasta.out)
		make_tree(make_msa.out.combine(make_dates_file.out))

	emit:
		align = make_msa.out
		tree = make_tree.out.tree

}

workflow create_ptree {
	take:
		ch_consensus_fasta
		ch_gene_database_fasta

	main:
		extract_sample_genes(ch_consensus_fasta.combine(ch_gene_database_fasta))
		make_multifasta(ch_gene_database_fasta.mix(extract_sample_genes.out).collect())
		make_dates_file(make_multifasta.out)
		make_msa(make_multifasta.out)
		make_tree(make_msa.out.combine(make_dates_file.out))

	emit:
		align = make_msa.out
		tree = make_tree.out.tree

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
	ch_blastdb_gtype_fasta = Channel.from(params.gtype_database)
	ch_blastdb_ptype_fasta = Channel.from(params.ptype_database)


	ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }

	main:
		//hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))

		// UNION DATABASE
		make_union_database(ch_blastdb_gtype_fasta, ch_blastdb_ptype_fasta)

		// QUALITY CONTROL 
		cutadapt(ch_fastq_input)
		fastp(cutadapt.out.trimmed_reads)
		fastp_json_to_csv(fastp.out.json)
		
		fastQC(fastp.out.trimmed_reads)
		fastq_check(fastp.out.trimmed_reads)
		fastq_check.out.formatted.collectFile(name: "${params.outdir}/fastq_qual/fastq_stats.tsv", keepHeader: true, skip: 1)
		
		// KRAKEN FILTERING
		run_kraken(fastp.out.trimmed_reads)
		kraken_filter(fastp.out.trimmed_reads.join(run_kraken.out.main))
		
		// DEHOSTING
		build_composite_reference(ch_human_ref.combine(make_union_database.out.fasta))
		index_composite_reference(build_composite_reference.out)
		dehost_fastq(
			kraken_filter.out.fastq,
			make_union_database.out.headers.first(), 
			index_composite_reference.out.first()
		) 

		// ASSEMBLY
		assembly(dehost_fastq.out.fastq)

		// ASSEMBLY QC
		run_quast(assembly.out.map{ it -> it[1]}.collect())

		// GENOTYPING / PTYPING 
		genotyping(assembly.out, ch_blastdb_gtype_fasta)
		ptyping(assembly.out, ch_blastdb_ptype_fasta)

		// Collect typing results 
		genotyping.out.filter.collectFile(name: "${params.outdir}/blastn/final/gtypes.tsv", keepHeader: true, skip: 1)
		ptyping.out.filter.collectFile(name: "${params.outdir}/blastn/final/ptypes.tsv", keepHeader: true, skip: 1)


		// Use large database to find reference
		if (params.reference_database){
			// SEARCH REFERENCE DB 
			find_reference(assembly.out)

			// SINGLE REFERENCE MAPPING 
			create_bwa_index(find_reference.out.ref)
			map_reads(dehost_fastq.out.fastq.join(create_bwa_index.out))
			sort_filter_index_sam(map_reads.out)

			ch_best_reference = find_reference.out.ref
			ch_bamfile = sort_filter_index_sam.out

		// Synthesize new reference from scratch
		} else {
			// COMBINE REFERENCES
			combine_references(genotyping.out.main.join(ptyping.out.main))
			combine_references.out.blast.collectFile(name: "${params.outdir}/blastn/final/final_types.tsv", keepHeader: true, skip: 1, newLine: true)

			// COMPETITIVE MAPPING 
			create_bwa_index(combine_references.out.ref)
			map_reads(dehost_fastq.out.fastq.join(create_bwa_index.out))
			sort_filter_index_sam(map_reads.out)

			// GENERATE NEW MERGED REFERENCE
			merge_fasta_bam(combine_references.out.ref.join(sort_filter_index_sam.out))

			ch_best_reference = merge_fasta_bam.out.ref
			ch_bamfile = merge_fasta_bam.out.bam
		}

		make_pileup(ch_best_reference.join(ch_bamfile))
		//create_pileup.out.metrics.collectFile(name: "${params.outdir}/qc/pileups", keepHeader: true, skip: 1)
		run_qualimap(ch_bamfile)

		// PLOT COVERAGE
		get_coverage(ch_bamfile)
		plot_coverage(get_coverage.out.coverage_file)

		// VARIANT CALLING
		create_fasta_index(ch_best_reference)
		run_freebayes(ch_bamfile.join(create_fasta_index.out))
		run_mpileup(ch_bamfile.join(create_fasta_index.out))
		get_common_snps(run_freebayes.out.join(run_mpileup.out))
		
		// CONSENSUS GENERATION 
		mask_low_coverage(ch_bamfile)
		make_consensus(get_common_snps.out.join(ch_best_reference).join(mask_low_coverage.out))

		create_gtree(make_consensus.out, genotyping.out.gene_db)
		create_ptree(make_consensus.out, ptyping.out.gene_db)

		run_custom_qc(ch_bamfile.join(ch_best_reference).join(make_consensus.out))
		run_custom_qc.out.csv.collectFile(name: "${params.outdir}/qc/custom/qc_all.csv", keepHeader: true, skip: 1)

		make_multifasta(make_consensus.out.map{it -> it[1]}.collect())
		make_msa(make_multifasta.out)
		make_tree(make_msa.out)

		// Collect all the relevant files for MULTIQC
		ch_multiqc_inputs = Channel.empty()
		ch_multiqc_inputs = ch_multiqc_inputs.mix(fastQC.out.zip.map{ it -> [it[1], it[2]]}.collect())
		ch_multiqc_inputs = ch_multiqc_inputs.mix(cutadapt.out.log)
		ch_multiqc_inputs = ch_multiqc_inputs.mix(run_kraken.out.report)
		ch_multiqc_inputs = ch_multiqc_inputs.mix(run_quast.out.tsv)
		ch_multiqc_inputs = ch_multiqc_inputs.mix(fastp.out.json.map{it -> it[1]})
		ch_multiqc_inputs = ch_multiqc_inputs.mix(run_qualimap.out.main)
		multiqc(ch_multiqc_inputs.collect().ifEmpty([]) )

		// ch_provenance = FluViewer.out.provenance
		// ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], [it[1]] << it[2]] }
		// ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(cutadapt.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
		// ch_provenance = ch_provenance.join(ch_fastq_input.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
		// collect_provenance(ch_provenance)

}