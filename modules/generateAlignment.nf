#!/bin/bash nextflow
//params.outdir = 'results'

include set_templates_path from './functions.nf'
path_templates = set_templates_path()

process REG_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method
    each bucket_size

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/regressive_align/reg_${align_method}.sh"   
}

process PROG_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.*.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/progressive_align/prog_${align_method}.sh"
}

process SLAVE_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method - $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method
    each bucket_size
    each slave_method   

    output:
    val align_method, emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.${align_method}.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/slave_align/slave_${align_method}.sh"
}

process DYNAMIC_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    val align_method
    each bucket_size
    each dynamicX
    path (dynamicConfig)
    tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize)
    val dynamicValues

    output:
    val "${dynamicValues}_${params.db}", emit: alignMethod
    val tree_method, emit: treeMethod
    val "${bucket_size}_${dynamicX}", emit: bucketSize
    tuple val (id), path("*.with.${tree_method}.tree.aln"), emit: alignmentFile 
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    // template can be default or using a config file
    template "${path_templates}/dynamic_align/dynamic_${align_method}.sh"
}

process POOL_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method
    each bucket_size

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.pool_${bucket_size}.${align_method}.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/pool_align/pool_${align_method}.sh"
}

process TCOFFEE_ALIGNER{
container 'fe94905825c8'
    tag "$tc_mode  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), path(seqs)
    val tc_mode

    output:
    val tc_mode, emit: tcMode
    tuple val (id), path ("${id}.tcoffee.*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:    
    if( tc_mode == 'default' )
        """
        t_coffee -seq $seqs > ${id}.tcoffee.${tc_mode}.aln
        """ 
    else if( tc_mode == 'quickaln' )
        """
        t_coffee -seq $seqs -mode quickaln > ${id}.tcoffee.${tc_mode}.aln
        """    
    else if( tc_mode == 'mcoffee' )
        """
        t_coffee -seq $seqs -mode mcoffee > ${id}.tcoffee.${tc_mode}.aln
        """ 
    else if( tc_mode == 'accurate' )
        """
        t_coffee -seq $seqs -mode accurate > ${id}.tcoffee.${tc_mode}.aln
        """            
    else if( tc_mode == 'fmcoffee' )
        """
        t_coffee -seq $seqs -mode fmcoffee > ${id}.tcoffee.${tc_mode}.aln
        """
    else if( tc_mode == 'psicoffee' )
        """
        t_coffee -seq $seqs -mode psicoffee > ${id}.tcoffee.${tc_mode}.aln
        """
    else if( tc_mode == 'expresso' )
        """
        t_coffee -seq $seqs -mode expresso -pdb_type dn > ${id}.tcoffee.${tc_mode}.aln
        """ 
    else if( tc_mode == 'procoffee' )
        """
        t_coffee -seq $seqs -mode procoffee > ${id}.tcoffee.${tc_mode}.aln
        """        
    else if( tc_mode == '3dcoffee' )
        """
        t_coffee -seq $seqs -method sap_pair -template_file sh3.template_file > ${id}.tcoffee.${tc_mode}.aln
        """
    else if( tc_mode == 'trmsd' )
        """
        t_coffee -seq $seqs -method mustang_pair -template_file crd.template_file

        t_coffee -other_pg trmsd -aln crd.aln -template_file crd.template_file
        """
    else if( tc_mode == 'rcoffee' )
        """
        t_coffee -seq $seqs -mode rcoffee -outfile RNA_rcoffee.aln > ${id}.tcoffee.${tc_mode}.aln
        """
    else
        error "Invalid alignment mode: ${tc_mode}"
}