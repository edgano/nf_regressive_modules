#!/bin/bash nextflow
params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process TREE_GENERATION {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method on $id"
    publishDir "${params.outdir}/trees", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(seqs)
    each tree_method

    output:
    tuple val (id), val(tree_method), path ("${id}.${tree_method}.dnd"), emit: trees

    script:
    template "${path_templates}/tree/tree_${tree_method}.sh"
}