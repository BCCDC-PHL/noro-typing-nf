process make_multifasta {
    publishDir "${params.outdir}/phylo/${custom_dir}/align", pattern: "*fasta" , mode:'copy'

    input: 
	path(sequences)

    output:
    path("${params.run_name}.multi.fasta")

    script:
    custom_dir = task.ext.custom_dir ?: 'full_genome'
    
	"""
	cat ${sequences} > temp.fasta
    filter_fasta.py temp.fasta ${params.run_name}.multi.fasta
    rm temp.fasta
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

process extract_sample_genes { 
    publishDir "${params.outdir}/phylo/${custom_dir}", pattern: "*.fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(consensus), path(blast_gene_db)

    output:
    path("${sample_id}*.fasta")

    script:
    gene = task.ext.gene ?: ''
    custom_dir = task.ext.custom_dir ?: 'full_genome'

	"""
    extract_genes.py sample --ref ${blast_gene_db} --query ${consensus} --gene ${gene} --output ${sample_id}_${gene}.fasta 
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

process make_annotated_tree {

    conda "${projectDir}/environments/plot.yaml"

    input:
    tuple path(tree)

    output:
    path("${params.run_name}_annotated_tree.pdf")

    script:

    """
    build_context_tree.R ${tree} ${params.run_name}_annotated_tree.pdf
    """
}