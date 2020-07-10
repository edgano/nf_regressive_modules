#!/bin/bash nextflow

include set_templates_path from './functions.nf'
path_templates = set_templates_path()

process SLAVE_CLUSTALO {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "clustalo", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.clustalo.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    for $slave in ${slave_method}
    do
        t_coffee -reg -reg_method clustalo_msa \
            -seq ${seqs} \
            -reg_tree ${guide_tree} \
            -child_tree $slave \
            -reg_nseq ${bucket_size} \
            -reg_homoplasy \
            -outfile ${id}.slave_${bucket_size}.clustalo.with.${tree_method}_${slave_method}.tree.aln   
    done
    """
}

process SLAVE_FAMSA{
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "famsa", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.famsa.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method famsa_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.famsa.with.${tree_method}_${slave_method}.tree.aln  
    """
}

process SLAVE_MAFFT_FFTNS1 {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "mafft-fftns1", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.mafft-fftns1.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method mafftfftns1_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.mafft-fftns1.with.${tree_method}_${slave_method}.tree.aln  
    """
}

process SLAVE_MAFFT_GINSI {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "mafft-ginsi", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.mafft-ginsi.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method mafftginsi_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.mafft-ginsi.with.${tree_method}_${slave_method}.tree.aln  
    """
}

process SLAVE_MAFFT_SPARSECORE {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "mafft-sparsecore", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.mafft-sparsecore.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method mafftsparsecore_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.mafft-sparsecore.with.${tree_method}_${slave_method}.tree.aln  
    """
}

process SLAVE_MAFFT {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "mafft", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.mafft.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method mafft_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.mafft.with.${tree_method}_${slave_method}.tree.aln  
    """
}

process SLAVE_MSAPROBS {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "msaProbs", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.msaProbs.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method msaprobs_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.msaProbs.with.${tree_method}_${slave_method}.tree.aln  
    """
}

process SLAVE_MUSCLE {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "muscle", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.muscle.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method muscle_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.muscle.with.${tree_method}_${slave_method}.tree.aln  
    """
}

process SLAVE_PROBCONS {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "probcons", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.probcons.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method probcons_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.probcons.with.${tree_method}_${slave_method}.tree.aln  
    """
}

process SLAVE_TCOFFEE {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "tcoffee", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.tcoffee.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method tcoffee_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.tcoffee.with.${tree_method}_${slave_method}.tree.aln  
    """
}

process SLAVE_UPP {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method & $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size
    each slave_method   

    output:
    val "upp", emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.upp.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method upp_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -child_tree ${slave_method} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.slave_${bucket_size}.upp.with.${tree_method}_${slave_method}.tree.aln  
    """
}