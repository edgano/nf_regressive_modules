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
    val (tree_method)
    //each tree_method from params.tree_method.tokenize(',') 

    output:
    tuple val(id), val (tree_method), file("${id}.${tree_method}.dnd") 

    script:
    template "${path_templates}/tree/tree_${tree_method}.sh"
}