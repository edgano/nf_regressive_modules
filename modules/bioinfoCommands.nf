#!/bin/bash nextflow

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process BLASTP {
    container 'edgano/blast_pdb'
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
    container 'edgano/blast_pdb'
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

process CONVERT2GAP{
    container '14c27fbd79b9'        // base with C
    tag "fasta2gap: ${fasta2gap} - $id"
    publishDir "${params.outdir}/aln2gap", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(infile)
    val(fasta2gap)                  //if true  : convert fasta to gap
                                    //if false : convert gap to fasta
    output:
    path "${id}_{F2G.gap,G2F.fa}"

    script:
    """
    ${baseDir}/bin/bioinfoCommands/gap ${infile} ${fasta2gap}

    if (${fasta2gap}); then               
        mv resultF2G.gap ${id}_F2G.gap
    else                              
        mv resultG2F.fa ${id}_G2F.fa
    fi
    """    
}