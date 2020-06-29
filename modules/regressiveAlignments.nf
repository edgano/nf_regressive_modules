#!/bin/bash nextflow
//params.outdir = 'results'

include set_templates_path from './functions.nf'
path_templates = set_templates_path()

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
         -method sap_pair mustang_pair \
         $template_filter $libs_filter \
         -outfile ${id}.reg_${bucket_size}.3dMcoffee.with.${tree_method}.tree.aln
    """
}