process make_multifasta {
    publishDir "${params.outdir}/phylo/align", pattern: "*fasta" , mode:'copy'

    input: 
	path(sequences)

    output:
    path("${params.run_name}.multi.fasta")

	"""
	cat ${sequences} > ${params.run_name}.multi.fasta
	"""

}

process make_msa {

    label 'ultra'

    publishDir "${params.outdir}/phylo/${custom_dir}/align", pattern: "*fasta" , mode:'copy'

    input: 
    path(multifasta)

    output:
    path("${params.run_name}.align.fasta")

    script: 
    custom_dir = task.ext.custom_dir ?: 'full_genome'
	"""
	mafft --auto --thread ${task.cpus} ${multifasta} > ${params.run_name}.align.fasta
	"""
}



process extract_genes_samples { 
    publishDir "${params.outdir}/phylo/polymerase", pattern: "*rdrp.fasta" , mode:'copy'
    publishDir "${params.outdir}/phylo/capsid", pattern: "*vp1.fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(consensus)

    output:
    path("${sample_id}*vp1.fasta"), emit: gtype
    path("${sample_id}*rdrp.fasta"), emit: ptype

	"""
    extract_genes.py --genes ${params.gene_positions} --gtype --ptype --prefix ${sample_id} ${consensus} ${params.g1_reference} 
	"""
}

process extract_genes_refs { 
    publishDir "${params.outdir}/phylo/${custom_dir}/", pattern: "*fasta" , mode:'copy'

    input: 
    path(reference)

    output:
    path("*.fasta")

    script: 
    // workflow_type = "${reference}" =~ /gtype/ ? "gtype" : "ptype" 
    custom_dir = task.ext.custom_dir ?: 'full_genome'
    workflow = task.ext.workflow ? "--${task.ext.workflow}" : ''
	"""
    extract_genes.py --genes ${params.gene_positions} ${workflow} --prefix ${reference.simpleName} ${reference} ${params.g1_reference} 
	"""
}

process make_tree {

    label 'ultra'

    publishDir "${params.outdir}/phylo/${custom_dir}/tree", pattern: "*nwk" , mode:'copy'

    input: 
    path(alignment)

    output:
    path("*nwk")

    script: 
    custom_dir = task.ext.custom_dir ?: 'full_genome'
    workflow = task.ext.workflow ? "-${task.ext.workflow}" : ''
    
	"""
    iqtree -T ${task.cpus} -m GTR -s ${alignment} --prefix ${params.run_name}${workflow}
    mkdir -p archive
    mv ${params.run_name}* archive
    mv archive/${params.run_name}${workflow}.treefile ./${params.run_name}${workflow}.nwk
	"""
}