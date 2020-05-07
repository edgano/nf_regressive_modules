#!/bin/bash nextflow
params.outdir = 'results'

process TREE_GENERATION {
    container 'edgano/homoplasy:latest'
    tag "$tree_method on $id"
    publishDir "${params.outdir}/trees"

    input:
    tuple id, path(seqs)
    val (tree_method)

    //each tree_method from tree_methods.tokenize(',') 

    output:
    file("${id}.${tree_method}.dnd") 

    script:
    template "${baseDir}/modules/regressive_alignment/templates/tree/tree_${tree_method}.sh"
}