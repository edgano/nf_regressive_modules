#!/bin/bash nextflow

include set_templates_path from './functions.nf'
path_templates = set_templates_path()

process BLASTP {
    container 'edgano/tcoffee:pdb'
    label 'process_high'
    tag "$id - $params.db"
    publishDir "${params.outdir}"

    input:
    tuple val(id), path(seqs)
    val (numThreads)

    output:
    //tuple val (id), path ("${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln"), emit: alignmentFile
    //path "${id}.homoplasy", emit: homoplasyFile
    //path ".command.trace", emit: metricFile

    script:
    """
    blastp -query ${seqs} -db ${params.database_path} -num_threads ${numThreads}
    """
}

process MAKEBLASTDB{
    container 'edgano/tcoffee:pdb'
    tag "$id - $params.db"
    publishDir "${params.outdir}/db", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(infile)
    val (inputType)
    val (dbType)

    output:
    //tuple val (id), path ("${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln"), emit: alignmentFile
    //path "${id}.homoplasy", emit: homoplasyFile
    //path ".command.trace", emit: metricFile

    script:
    """
    makeblastdb -in ${infile} -input_type ${inputType} -dbtype ${dbType} 
    """
}