#!/bin/bash nextflow

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process PROG_CLUSTALO {
    container 'edgano/tcoffee:pdb'
    tag "clustalo - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "clustalo", emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.clustalo.with.${tree_method}.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    clustalo --infile=${seqs} \
        --guidetree-in=${guide_tree} \
        --outfmt=fa \
        -o ${id}.prog.clustalo.with.${tree_method}.tree.aln 
    """
}

process PROG_FAMSA {
    container 'edgano/tcoffee:pdb'
    tag "famsa - $tree_method  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "famsa", emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.famsa.with.${tree_method}.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    famsa -gt import ${guide_tree} ${seqs} ${id}.prog.famsa.with.${tree_method}.tree.aln
    """
}

process PROG_MAFFT_FFTNS1 {
    container 'edgano/tcoffee:pdb'
    tag "mafft-fftns1 - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "mafft-fftns1", emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.mafft-fftns1.with.${tree_method}.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1

    t_coffee -other_pg seq_reformat -in ${guide_tree} -input newick -in2 ${seqs} -input2 fasta_seq -action +newick2mafftnewick >> ${id}.mafftnewick
    newick2mafft.rb 1.0 ${id}.mafftnewick > ${id}.mafftbinary
    mafft --retree 1 --anysymbol --treein ${id}.mafftbinary ${seqs} > ${id}.prog.mafft-fftns1.with.${tree_method}.tree.aln
    """
}

process PROG_MAFFT_GINSI {
    container 'edgano/tcoffee:pdb'
    tag "mafft-ginsi - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "mafft-ginsi", emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.mafft-ginsi.with.${tree_method}.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -other_pg seq_reformat -in ${guide_tree} -input newick -in2 ${seqs} -input2 fasta_seq -action +newick2mafftnewick >> ${id}.mafftnewick
    newick2mafft.rb 1.0 ${id}.mafftnewick > ${id}.mafftbinary
    ginsi --treein ${id}.mafftbinary ${seqs} > ${id}.prog.mafft-ginsi.with.${tree_method}.tree.aln
    """
}

process PROG_MAFFT_SPARSECORE {
    container 'edgano/tcoffee:pdb'
    tag "mafft-sparsecore - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "mafft-sparsecore", emit: alignMethod
    val "NO_TREE", emit: treeMethod
    tuple val (id), path ("${id}.prog.mafft-sparsecore.with.NO_TREE.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    replace_U.pl ${seqs}

    mafft-sparsecore.rb -i ${seqs} > ${id}.prog.mafft-sparsecore.with.NO_TREE.tree.aln
    """
}

process PROG_MAFFT {
    container 'edgano/tcoffee:pdb'
    tag "mafft - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "mafft", emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.mafft.with.${tree_method}.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -other_pg seq_reformat -in ${guide_tree} -input newick -in2 ${seqs} -input2 fasta_seq -action +newick2mafftnewick >> ${id}.mafftnewick

    newick2mafft.rb 1.0 ${id}.mafftnewick > ${id}.mafftbinary

    mafft --anysymbol --treein ${id}.mafftbinary ${seqs} > ${id}.prog.mafft.with.${tree_method}.tree.aln
    """
}

process PROG_MSAPROBS {
    container 'edgano/tcoffee:pdb'
    tag "msaprobs - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "msaprobs", emit: alignMethod
    val "NO_TREE", emit: treeMethod
    tuple val (id), path ("${id}.prog.msaprobs.with.NO_TREE.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    msaprobs ${seqs} -o ${id}.prog.msaprobs.with.NO_TREE.tree.aln
    """
}

process PROG_MUSCLE {
    container 'edgano/tcoffee:pdb'
    tag "muscle - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "muscle", emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.muscle.with.${tree_method}.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    muscle -in ${seqs} -out ${id}.prog.muscle.with.${tree_method}.tree.aln -usetree_nowarn ${guide_tree}
    """
}

process PROG_PROBCONS {
    container 'edgano/tcoffee:pdb'
    tag "probcons - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "probcons", emit: alignMethod
    val "NO_TREE", emit: treeMethod
    tuple val (id), path ("${id}.prog.probcons.with.NO_TREE.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    probcons ${seqs} > ${id}.prog.probcons.with.NO_TREE.tree.aln 
    """
}

/*process PROG_TCOFFEE {
    container 'edgano/tcoffee:pdb'
    tag "tcoffee - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "tcoffee", emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.tcoffee.with.${tree_method}.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -seq ${seqs} -tree ${guide_tree} \
         -outfile ${id}.prog.tcoffee.with.${tree_method}.tree.aln 
    """
}*/

process PROG_UPP {
    container 'edgano/tcoffee:pdb'
    tag "upp - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)

    output:
    val "upp", emit: alignMethod
    val "NO_TREE", emit: treeMethod
    tuple val (id), path ("${id}.prog.upp.with.NO_TREE.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    replace_U.pl ${seqs}

    run_upp.py -s ${seqs} \
            -m amino \
            -x 1 \
            -o ${id}.prog.upp

    mv ${id}.prog.upp_alignment.fasta ${id}.prog.upp.with.NO_TREE.tree.aln
    """
}