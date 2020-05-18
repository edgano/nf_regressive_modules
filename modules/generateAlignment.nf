#!/bin/bash nextflow
//params.outdir = 'results'

moduleDir="$baseDir/modules/"
path_templates = "${moduleDir}/templates"

process REG_ALIGNER {
    container 'edgano/tcoffee:psi'
    tag "$align_method - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments"

    input:
    tuple id, path(seqs)
    each align_method
    each bucket_size
    tuple val(id), val(tree_method), file(guide_tree)

    output:
    val id, emit:id
    path "${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln", emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/regressive_align/reg_${align_method}.sh"
}

process PROG_ALIGNER {
    container 'edgano/tcoffee:psi'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments"

    input:
    tuple id, path(seqs)
    each (align_method)
    tuple val(id), val(tree_method), file(guide_tree)

    output:
    val id, emit:id
    path "${id}.prog.${align_method}.with.${tree_method}.tree.aln", emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/progressive_align/prog_${align_method}.sh"
}

process SLAVE_ALIGNER {
    container 'edgano/tcoffee:psi'
    tag "$align_method - $tree_method - $slave_method on $id"
    publishDir "${params.outdir}/alignments"

    input:
    tuple id, path(seqs)
    each (align_method)
    each bucket_size
    tuple val(id), val(tree_method), file(guide_tree)
    val (slave_method)

    output:
    val id, emit:id
    path "${id}.slave_${bucket_size}.${align_method}.with.${tree_method}_${slave_method}.tree.aln", emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile
    val "${tree_method}_${slave_method}", emit:tree_method

    script:
    template "${path_templates}/slave_align/slave_${align_method}.sh"
}

process DYNAMIC_ALIGNER {
    container 'edgano/tcoffee:psi'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments"

    input:
    tuple id, path(seqs)
    each (align_method)
    each bucket_size
    tuple val(id), val(tree_method), file(guide_tree)

    output:
    file("${id}.dynamic_${bucket_size}.dynamicSize.${dynamic_size}.${align_method}.with.${tree_method}.tree.aln") 
    path ".command.trace", emit: metricFile

    script:
    // template "${path_templates}/dynamic_align/dynamic_${align_method}.sh"
    // the above template is not declared yet, thus I call the following one
    template "${path_templates}/dynamic_align/dynamic_DEFAULT.sh"
}
process POOL_ALIGNER {
    container 'edgano/tcoffee:psi'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments"

    input:
    tuple id, path(seqs)
    each (align_method)
    each bucket_size
    tuple val(id), val(tree_method), file(guide_tree)

    output:
    val id, emit:id
    path "${id}.pool_${bucket_size}.${align_method}.with.${tree_method}.tree.aln", emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/pool_align/pool_${align_method}.sh"
}