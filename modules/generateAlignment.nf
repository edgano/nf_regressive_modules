#!/bin/bash nextflow
//params.outdir = 'results'

moduleDir="$baseDir/modules/"
path_templates = "${moduleDir}/templates"

process REG_ALIGNER {
    container 'edgano/tcoffee:psi'
    tag "$align_method - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple id, path(seqs)
    each align_method
    each bucket_size
    each tree_method
    path (guide_tree)

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    //TODO refactor emit alignmentFile

    script:
    template "${path_templates}/regressive_align/reg_${align_method}.sh"
}

process PROG_ALIGNER {
    container 'edgano/tcoffee:psi'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple id, path(seqs)
    each align_method
    each tree_method
    path (guide_tree)

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.${align_method}.with.${tree_method}.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/progressive_align/prog_${align_method}.sh"
}

process SLAVE_ALIGNER {
    container 'edgano/tcoffee:psi'
    tag "$align_method - $tree_method - $slave_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple id, path(seqs)
    each align_method
    each bucket_size
    each tree_method
    path guide_tree
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
    tuple id, path(seqs)
    each align_method
    each bucket_size
    each dynamic_size
    each tree_method
    path (guide_tree)

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    val "${bucket_size}_${dynamic_size}", emit: bucketSize
    tuple val (id), path("${id}.dynamic_${bucket_size}.dynamicSize.${dynamic_size}.${align_method}.with.${tree_method}.tree.aln"), emit: alignmentFile 
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    // template "${path_templates}/dynamic_align/dynamic_${align_method}.sh"
    // the above template is not declared yet, thus I call the following one
    template "${path_templates}/dynamic_align/dynamic_DEFAULT.sh"
}

process POOL_ALIGNER {
    container 'edgano/tcoffee:psi'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple id, path(seqs)
    each align_method
    each bucket_size
    each tree_method
    path (guide_tree)

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