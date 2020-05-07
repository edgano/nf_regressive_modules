#!/bin/bash nextflow
params.outdir = 'results'

process COMBINE_SEQS {
    container 'edgano/homoplasy:latest'
    tag "COMBINE SEQ on $id"
    publishDir params.outdir

    input:
    tuple id, path(seqs)
    tuple id2, path(refs)

    output:
    file("${id}.fa") 

    script:
    """
    # CREATE A FASTA FILE CONTAINING ALL SEQUENCES (SEQS + REFS)
    t_coffee -other_pg seq_reformat -in ${refs} -output fasta_seq -out refs.tmp.fa
    t_coffee -other_pg seq_reformat -in ${seqs} -output fasta_seq -out seqs.tmp.fa
    cat refs.tmp.fa > completeSeqs.fa
    cat seqs.tmp.fa >> completeSeqs.fa
    # SHUFFLE ORDER OF SEQUENCES
    t_coffee -other_pg seq_reformat -in completeSeqs.fa -output fasta_seq -out ${id}.fa -action +reorder random
    """
}