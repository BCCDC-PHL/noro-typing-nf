process make_multifasta {
    publishDir "${params.outdir}/phylo/${custom_dir}/sequences", pattern: "*fasta" , mode:'copy'

    input: 
	path(sequences)

    output:
    path("${params.run_name}_${workflow}_multi.fasta")

    script:
    custom_dir = task.ext.custom_dir ?: 'full_genome'
    workflow = task.ext.workflow ?: ''
    
	"""
	cat ${sequences} > temp.fasta
    filter_fasta.py main temp.fasta ${params.run_name}_${workflow}_multi.fasta
    rm temp.fasta
	"""

}
process get_background_sequences {
    publishDir "${params.outdir}/phylo/${custom_dir}/sequences", pattern: "*fasta" , mode:'copy'

    input: 
	path(multifasta)
    path(results_path)
    
    output:
    path("${params.run_name}_${custom_dir}_bg.fasta")

    script:
    custom_dir = task.ext.custom_dir ?: 'full'
    workflow = task.ext.workflow ?: 'NONE'
    
	"""
    get_background_seqs.py --infasta ${multifasta} --gene ${custom_dir} --outfasta ${params.run_name}_${custom_dir}_bg.fasta ${results_path}
	"""
}

process make_msa {

    label 'ultra'

    publishDir "${params.outdir}/phylo/${custom_dir}/align", pattern: "*fasta" , mode:'copy'

    input: 
    path(fastas)

    output:
    path("${params.run_name}.align.fasta")

    script: 
    custom_dir = task.ext.custom_dir ?: 'full_genome'
	"""
    printf -- "- process_name: mafft\\n" > mafft_provenance.yml
    printf -- "  tool_name: mafft\\n  tool_version: \$(mafft --version 2>&1 )\\n" >> mafft_provenance.yml

    cat ${fastas} > multi.fasta &&
	mafft --auto --thread ${task.cpus} multi.fasta > ${params.run_name}.align.fasta
	"""
}

process extract_sample_genes { 

	tag { sample_id }
    
    errorStrategy 'ignore'

    publishDir "${params.outdir}/phylo/${custom_dir}", pattern: "*.fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(consensus), path(blast_gene_db), path(blast_gene_pos)

    output:
    path("${sample_id}*.fasta")

    script:
    gene = task.ext.gene ?: ''
    custom_dir = task.ext.custom_dir ?: 'full_genome'

	"""
    extract_genes.py sample --ref ${blast_gene_db} --query ${consensus} --positions ${blast_gene_pos} --gene ${gene} --outfasta ${sample_id}_${gene}.fasta 
	"""
}

process make_dates_file { 
    publishDir "${params.outdir}/phylo/${custom_dir}/tree", pattern: "*" , mode:'copy'

    input: 
    path(fastas)

    output:
    path("${out_name}")

    script:
    gene = task.ext.gene ?: ''
    workflow = task.ext.workflow ?: 'global'
    custom_dir = task.ext.custom_dir ?: 'global'
    out_name = "${params.run_name}_${workflow}_dates.tsv"
	"""
    cat ${fastas} > multi.fasta &&
    get_dates.py multi.fasta --output ${out_name}
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
    printf -- "- process_name: iqtree\\n" > iqtree_provenance.yml
    printf -- "  tool_name: iqtree\\n  tool_version: \$(iqtree --version 2>&1 | head -n1 | cut -d' ' -f4)\\n" >> iqtree_provenance.yml


    iqtree -T ${task.cpus} -m GTR -s ${alignment} ${dates} --prefix ${params.run_name}_${workflow}
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