#!/bin/bash nextflow
params.outdir = 'results'

process REG_ALIGNER {
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments"

    input:
    tuple id, path(seqs)
    val (align_method)
    val (tree_method)
    val (bucket_size)
    file (guide_tree)
//  each bucket_size from params.buckets.tokenize(',')

    output:
    val id, emit:id
    path "${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln", emit: alignment

    script:
    template "regressive_align/reg_${align_method}.sh"
}

process PROG_ALIGNER {
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments"

    input:
    tuple id, path(seqs)
    val (align_method)
    val (tree_method)
    file (guide_tree)

    output:
    file("${id}.prog.${align_method}.with.${tree_method}.tree.aln") 

    script:
    template "${baseDir}/modules/regressive_alignment/templates/progressive_align/prog_${align_method}.sh"
}

process SLAVE_ALIGNER {
    tag "$align_method - $tree_method - $slave_method on $id"
    publishDir "${params.outdir}/alignments"

    input:
    tuple id, path(seqs)
    val (align_method)
    val (tree_method)
    val (slave_method)
    file (guide_tree)

    output:
    file("${id}.slave.${align_method}.with.${tree_method}.tree.slave.${slave_method}.aln") 

    script:
    template "${baseDir}/modules/regressive_alignment/templates/slave_align/slave_${align_method}.sh"
}

process DYNAMIC_ALIGNER {
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments"

    input:
    tuple id, path(seqs)
    val (align_method)
    val (tree_method)
    file (guide_tree)

    output:
    file("${id}.dynamic.${bucket_size}.dynamicSize.${dynamic_size}.${align_method}.with.${tree_method}.tree.aln") 

    script:
    template "${baseDir}/modules/regressive_alignment/templates/dynamic_align/dynamic_${align_method}.sh"
}