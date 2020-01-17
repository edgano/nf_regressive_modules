#!/bin/bash nextflow
params.outdir = 'results'

process COMBINE_SEQS {
    tag "COMBINE SEQ on $id"
    publishDir params.outdir

    input:
    tuple id, path(reads)

    output:
    file("${id}.fa") 

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}