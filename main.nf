#!/usr/bin/env nextflow

/*
== V0.1  ==
*/

import java.time.LocalDateTime

nextflow.enable.dsl = 2

// include { pipeline_provenance } from './modules/provenance.nf'
// include { collect_provenance } from './modules/provenance.nf'
include { fastQC; fastq_check; run_quast; run_qualimap; run_custom_qc} from './modules/qc.nf'
include { make_union_database; cutadapt; fastp; run_kraken } from './modules/prep.nf' 
include { build_composite_reference; dehost_fastq } from './modules/prep.nf'
include { prep_database; make_blast_database; extract_genes_blast; run_self_blast; run_blastn; run_blastx; combine_references} from './modules/blast.nf'
//include { p_make_blast_database; p_run_self_blast; p_run_blastn; p_filter_alignments; p_get_best_references; p_run_blastx } from './modules/p_blast.nf'
include { assembly } from './modules/assembly.nf'
include { create_fasta_index; map_reads; sort_filter_index_sam; merge_fasta_bam} from './modules/mapping.nf'
include { run_freebayes; run_mpileup ; get_common_snps } from './modules/variant_calling.nf'
include { get_coverage; plot_coverage; make_pileup} from "./modules/coverage.nf"
include { mask_low_coverage; make_consensus } from './modules/consensus.nf'
include { make_multifasta; get_background_sequences; make_msa; make_dates_file; make_tree; extract_sample_genes} from './modules/phylo.nf'
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

workflow create_gtree {
	take:
		ch_consensus_fasta
		ch_sample_db_pos
		ch_gene_db
	main:
		extract_sample_genes(ch_consensus_fasta.join(ch_sample_db_pos))
		make_multifasta(extract_sample_genes.out.mix(ch_gene_db).collect()).set{ch_sequences}

		if (params.results_path){
			get_background_sequences(ch_sequences, Channel.fromPath(params.results_path))
			ch_sequences = ch_sequences.mix(get_background_sequences.out).collect()
		} 

		make_dates_file(ch_sequences)
		make_msa(ch_sequences)
		make_tree(make_msa.out.combine(make_dates_file.out))

	emit:
		align = make_msa.out
		tree = make_tree.out.tree

}

workflow create_ptree {
	take:
		ch_consensus_fasta
		ch_sample_db_pos
		ch_gene_db

	main:
		extract_sample_genes(ch_consensus_fasta.join(ch_sample_db_pos))
		make_multifasta(extract_sample_genes.out.mix(ch_gene_db).collect()).set{ch_sequences}

		if (params.results_path){
			get_background_sequences(ch_sequences, Channel.fromPath(params.results_path))
			ch_sequences = ch_sequences.mix(get_background_sequences.out).collect()
		} 

		make_dates_file(ch_sequences)
		make_msa(ch_sequences)
		make_tree(make_msa.out.combine(make_dates_file.out))

	emit:
		align = make_msa.out
		tree = make_tree.out.tree

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

		if (params.blast_refsearch_metric == 'bsr'){
			error "ERROR BLAST Reference search cannot use bsr (blast score ratio)."
		}

		sample_list = ch_fastq_input.map{it -> it[0]}.collect().first()

		// QUALITY CONTROL 
		cutadapt(ch_fastq_input)
		fastp(cutadapt.out.trimmed_reads)
		
		fastQC(fastp.out.trimmed_reads)
		fastq_check(fastp.out.trimmed_reads)
		fastq_check.out.formatted.collectFile(name: "${params.outdir}/fastq_qual/fastq_stats.tsv", keepHeader: true, skip: 1)
		
		// KRAKEN FILTERING
		run_kraken(fastp.out.trimmed_reads)
		
		// DEHOSTING
		make_union_database(ch_blastdb_gtype_fasta, ch_blastdb_ptype_fasta).first()
		build_composite_reference(ch_human_ref.combine(make_union_database.out.fasta))
		dehost_fastq(
			run_kraken.out.fastq,
			make_union_database.out.headers.first(), 
			build_composite_reference.out.first()
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

		// Synthesize new reference from scratch
		// COMBINE REFERENCES
		combine_references(genotyping.out.main.join(ptyping.out.main))
		//combine_references.out.blast.collectFile(name: "${params.outdir}/blastn/final/final_types.tsv", keepHeader: true, skip: 1, newLine: true)

		// COMPETITIVE MAPPING 
		map_reads(dehost_fastq.out.fastq.join(combine_references.out))
		sort_filter_index_sam(map_reads.out)

		// GENERATE NEW MERGED REFERENCE
		merge_fasta_bam(combine_references.out.join(sort_filter_index_sam.out))

		// CALCULATE COVERAGE OF THE SYNTHETIC METHOD 
		get_coverage(merge_fasta_bam.out.bam)

		ch_synthetic = get_coverage.out.metric
					.join(merge_fasta_bam.out.ref)
					.join(merge_fasta_bam.out.bam)
					.join(get_coverage.out.main)
					.combine(genotyping.out.db_pos)
					.combine(ptyping.out.db_pos)

		// If parameter passed, use comprehensive reference database to find reference
		if (params.reference_db_full){

			// reformat the synthetic channel for comparison purposes 
			ch_synthetic = ch_synthetic.map{ it -> [it[0], it.subList(1,it.size())]}
			
			// perform a BLAST search on the global reference database 
			global_reference_search(assembly.out, dehost_fastq.out.fastq)

			// make a single comparison channel to compare outputs between both methods 
			ch_compare = global_reference_search.out.main.join(ch_synthetic, remainder: true)

			outfile = file("${params.outdir}/qc/choices.txt")

			ch_compare.map{ name, vals1, vals2 ->
				if (vals2 == null){
					outfile.append("${name},${0},${vals1[0]}\n")
					[name] + vals1.subList(1, vals1.size())
				}else if (vals1[0] >= vals2[0]){
					outfile.append("${name},${vals2[0]},${vals1[0]}\n")
					[name] + vals1.subList(1, vals1.size())
				}else {
					outfile.append("${name},${vals2[0]},${vals1[0]}\n")
					[name] + vals2.subList(1, vals2.size())
				}
			}.set{ch_best_coverage}

		}else {
			ch_best_coverage = ch_synthetic.map{it -> [it[0]] + it.subList(2,it.size())}
		}

		ch_best_reference = ch_best_coverage.map{it -> [it[0], it[1]]}
		ch_bamfile = ch_best_coverage.map{it -> [it[0], it[2], it[3]]}
		ch_coverage = ch_best_coverage.map{it -> [it[0], it[4]]}
		ch_gene_db_vp1 = ch_best_coverage.map{it -> [it[0], it[5], it[6]]}
		ch_gene_db_rdrp = ch_best_coverage.map{it -> [it[0], it[7], it[8]]}

		make_pileup(ch_best_reference.join(ch_bamfile))
		//create_pileup.out.metrics.collectFile(name: "${params.outdir}/qc/pileups", keepHeader: true, skip: 1)
		run_qualimap(ch_bamfile)

		// PLOT COVERAGE
		// get_coverage(ch_bamfile)
		plot_coverage(ch_coverage)

		// VARIANT CALLING
		create_fasta_index(ch_best_reference)
		run_freebayes(ch_bamfile.join(create_fasta_index.out))
		run_mpileup(ch_bamfile.join(create_fasta_index.out))
		get_common_snps(run_freebayes.out.join(run_mpileup.out))
		
		// CONSENSUS GENERATION 
		mask_low_coverage(ch_bamfile)
		make_consensus(get_common_snps.out.join(ch_best_reference).join(mask_low_coverage.out))

		// GENE-SPECIFIC PHYLOGENETIC TREES
		create_gtree(make_consensus.out, ch_gene_db_vp1, ch_blastdb_gtype_fasta)
		create_ptree(make_consensus.out, ch_gene_db_rdrp, ch_blastdb_ptype_fasta)

		// CUSTOM QC 
		run_custom_qc(ch_bamfile.join(ch_best_reference).join(make_consensus.out))
		run_custom_qc.out.csv.collectFile(name: "${params.outdir}/qc/custom/qc_all.csv", keepHeader: true, skip: 1, newLine:true)

		make_multifasta(make_consensus.out.map{it -> it[1]}.collect()).set{ch_sequences}

		if (params.results_path){
			get_background_sequences(ch_sequences, Channel.fromPath(params.results_path))
			ch_sequences = ch_sequences.mix(get_background_sequences.out).collect()
		} 

		make_msa(ch_sequences)
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