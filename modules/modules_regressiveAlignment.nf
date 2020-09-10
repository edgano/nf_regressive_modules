#!/bin/bash nextflow

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process REG_CLUSTALO {
    container 'edgano/tcoffee:pdb'
    tag "Regressive clustalo - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "clustalo", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.clustalo.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method clustalo_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.clustalo.with.${tree_method}.tree.aln
    """
}

process REG_FAMSA{
    container 'edgano/tcoffee:pdb'
    tag "Regressive famsa - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "famsa", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.famsa.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method famsa_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.famsa.with.${tree_method}.tree.aln
    """
}

process REG_MAFFT_FFTNS1 {
    container 'edgano/tcoffee:pdb'
    tag "Regressive fftns1 - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "mafft-fftns1", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.mafft-fftns1.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method mafftfftns1_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.mafft-fftns1.with.${tree_method}.tree.aln
    """
}

process REG_MAFFT_GINSI {
    container 'edgano/tcoffee:protocols'
    tag "Regressive ginsi - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "mafft-ginsi", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.mafft-ginsi.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method mafftginsi_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.mafft-ginsi.with.${tree_method}.tree.aln
    """
}

process REG_MAFFT_SPARSECORE {
    container 'edgano/tcoffee:protocols'
    tag "Regressive sparsecore - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "mafft-sparsecore", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.mafft-sparsecore.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method mafftsparsecore_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.mafft-sparsecore.with.${tree_method}.tree.aln
    """
}

process REG_MAFFT {
    container 'edgano/tcoffee:protocols'
    tag "Regressive mafft - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "mafft", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.mafft.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method mafft_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.mafft.with.${tree_method}.tree.aln
    """
}

process REG_MSAPROBS {
    container 'edgano/tcoffee:protocols'
    tag "Regressive msaProbs - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "msaProbs", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.msaProbs.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export NO_MAFFT_BINARIES=1
    export MAFFT_BINARIES=''

    t_coffee -reg -reg_method msaprobs_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.msaProbs.with.${tree_method}.tree.aln
    """
}

process REG_MUSCLE {
    container 'edgano/tcoffee:protocols'
    tag "Regressive muscle - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "muscle", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.muscle.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method muscle_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.muscle.with.${tree_method}.tree.aln
    """
}

process REG_PROBCONS {
    container 'edgano/tcoffee:protocols'
    tag "Regressive probcons - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "probcons", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.probcons.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method probcons_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.probcons.with.${tree_method}.tree.aln
    """
}

process REG_PSICOFFEE {
    container 'edgano/tcoffee:protocols'
    tag "Regressive psicoffee - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "psicoffee", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.psicoffee.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    export blast_server_4_CLTCOFFEE=LOCAL
    export protein_db_4_CLTCOFFEE=${params.database_path}
    export VERBOSE_4_DYNAMIC=1

    t_coffee -reg -reg_method psicoffee_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -method sap_pair  \
         -outfile ${id}.reg_${bucket_size}.psicoffee.with.${tree_method}.tree.aln
    """
}

process REG_TCOFFEE {
    container 'edgano/tcoffee:protocols'
    tag "Regressive tcoffee - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "tcoffee", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.tcoffee.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method tcoffee_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.tcoffee.with.${tree_method}.tree.aln
    """
}

process REG_UPP {
    container 'edgano/tcoffee:protocols'
    tag "Regressive upp - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each bucket_size

    output:
    val "upp", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.upp.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -reg -reg_method upp_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.upp.with.${tree_method}.tree.aln
    """
}

process REG_3DALIGN {
    container 'edgano/tcoffee:protocols'
    tag "Regressive 3d_align - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(guide_tree), path(seqs), file(template), file(library)
    each bucket_size

    output:
    val "3dcoffee", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.3dcoffee.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    //ERROR TMalign_pair
    def template_filter = template.name != 'input.4' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.5' ? "-lib $library" : ''
    """
    t_coffee -reg -reg_method 3dcoffee_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -method sap_pair  \
         $template_filter $libs_filter \
         -outfile ${id}.reg_${bucket_size}.3dcoffee.with.${tree_method}.tree.aln
    """ 
}

process REG_3DMALIGN {
    container 'edgano/tcoffee:protocols'
    tag "Regressive 3dM_align - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(guide_tree), path(seqs), file(template), file(library)
    each bucket_size

    output:
    val "3dMcoffee", emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.3dMcoffee.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.4' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.5' ? "-lib $library" : ''
    """
    t_coffee -reg -reg_method 3dcoffee_msa \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -method sap_pair  \
         $template_filter $libs_filter \
         -outfile ${id}.reg_${bucket_size}.3dMcoffee.with.${tree_method}.tree.aln
    """
}