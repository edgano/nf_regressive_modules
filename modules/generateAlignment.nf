#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process REG_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'
    //publishDir "${params.outdir}/templates", pattern: '*.template_list'
    //publishDir "${params.outdir}/templates", pattern: '*.prf'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree)
    each align_method
    each bucket_size

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.*.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile
    path "time.txt", emit: timeFile
    //path "*.template_list", emit: templateFile
    //path "*.prf", emit: templateProfile
    
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
    label 'process_medium'

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

process TCOFFEE_ALIGNER {
    container 'edgano/tcoffee:protocols'
    tag "$tc_mode  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'
    //publishDir "${params.cache_path}", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), file(template), file(library)
    each tc_mode         

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile
    path "*.aln" optional true          //TCOFFEE.out[2].view()


    script:
    template "${path_templates}/tcoffee_align/tcoffee_${tc_mode}.sh"
}
