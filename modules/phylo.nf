process make_multifasta {
    publishDir "${params.outdir}/phylo/${custom_dir}/align", pattern: "*fasta" , mode:'copy'

    input: 
	path(sequences)

    output:
    path("${params.run_name}_${workflow}.multi.fasta")

    script:
    custom_dir = task.ext.custom_dir ?: 'full_genome'
    workflow = task.ext.workflow ?: ''
    
	"""
	cat ${sequences} > temp.fasta
    filter_fasta.py temp.fasta ${params.run_name}_${workflow}.multi.fasta
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

	tag { sample_id }
    
    errorStrategy 'ignore'

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

process make_dates_file { 
    publishDir "${params.outdir}/phylo/${custom_dir}/tree", pattern: "*" , mode:'copy'

    input: 
    path(multifasta)

    output:
    path("${out_name}")

    script:
    gene = task.ext.gene ?: ''
    custom_dir = task.ext.custom_dir ?: 'full_genome'
    out_name = "${multifasta.simpleName}_dates.tsv"
	"""
    get_dates.py ${multifasta} --output ${out_name}
	"""
}

process make_tree {

    label 'ultra'

    publishDir "${params.outdir}/phylo/${custom_dir}/tree", pattern: "*nwk" , mode:'copy'

    input: 
    path(infiles)

    output:
    path("*nwk"), emit: tree
    path("${params.run_name}_iqtree"), emit: archive

    script: 
    custom_dir = task.ext.custom_dir ?: 'full_genome'
    workflow = task.ext.workflow ?: 'full'
    alignment = infiles.size() == 2 ? infiles[0] : infiles
    dates = infiles.size() == 2 ? "--date ${infiles[1]}" : ''
	"""
    iqtree -T ${task.cpus} -m GTR -s ${alignment} --date ${dates} --prefix ${params.run_name}_${workflow}
    mkdir -p ${params.run_name}_iqtree
    mv ${params.run_name}_${workflow}* ${params.run_name}_iqtree
    cp ${params.run_name}_iqtree/${params.run_name}_${workflow}.treefile ./${params.run_name}_${workflow}.nwk
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