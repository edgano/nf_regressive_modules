#!/bin/bash nextflow
params.outdir = 'results'

moduleDir="$baseDir/modules/"
path_templates = "${moduleDir}/templates"

process TREE_GENERATION {
    container 'edgano/homoplasy:latest'
    tag "$tree_method on $id"
    publishDir "${params.outdir}/trees"

    input:
    tuple id, path(seqs)
    each tree_method

    output:
    val id, emit:id
    val tree_method, emit:treeMethod
    path "${id}.${tree_method}.dnd", emit:guideTree

    script:
    template "${path_templates}/tree/tree_${tree_method}.sh"
}