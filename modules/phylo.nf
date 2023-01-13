process make_multifasta {
    publishDir "${params.outdir}/phylo/align", pattern: "*fasta" , mode:'copy'

    input: 
	path(consensus)

    output:
    path("${params.run_name}.multi.fasta")

	"""
	cat ${consensus} > ${params.run_name}.multi.fasta
	"""

}

process add_background_sequences {
    // publishDir "${params.outdir}/phylo/align", pattern: "*fasta" , mode:'copy'

    input: 
	tuple path(multifasta), path(consensus_path)

    output:
    path("${params.run_name}.multi.fasta")

	"""
	cat ${consensus} > ${params.run_name}.multi.fasta
	"""
}

process make_msa {

    label 'ultra'

    publishDir "${params.outdir}/phylo/align", pattern: "*fasta" , mode:'copy'

    input: 
    path(multifasta)

    output:
    path("${params.run_name}.align.fasta")

	"""
	mafft --auto --thread ${task.cpus} ${multifasta} > ${params.run_name}.align.fasta
	"""
}

process make_tree {

    label 'ultra'

    publishDir "${params.outdir}/phylo/tree", pattern: "*" , mode:'copy'

    input: 
    path(alignment)

    output:
    path("tree_out*")

	"""
    iqtree -T ${task.cpus} -m GTR -s ${alignment} --prefix tree_out
	"""
}
