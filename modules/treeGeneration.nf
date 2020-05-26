#!/bin/bash nextflow
params.outdir = 'results'

moduleDir="$baseDir/modules/"
path_templates = "${moduleDir}/templates"

process TREE_GENERATION {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method on $id"
    publishDir "${params.outdir}/trees"

    input:
    tuple id, path(seqs)
    each tree_method

    output:
    tuple val (id), val(tree_method), path ("${id}.${tree_method}.dnd")

    script:
    template "${path_templates}/tree/tree_${tree_method}.sh"
}