#!/bin/bash nextflow
params.outdir = 'results'

path_templates = projectDir.name == "nf-benchmark" ? "$baseDir/modules/regressive_alignment/templates" : "$baseDir/templates"

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
    template "${path_templates}/tree/tree_${tree_method}.sh"
}