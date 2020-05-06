#!/bin/bash nextflow
params.outdir = 'results'

process EVAL_ALIGNMENT {
    tag "EVAL_ALIGNMENT on $id"
    publishDir "${params.outdir}/score"

    input:
    file (test_alignment)
    tuple id, file(ref_alignment)
    val (align_method)
    val (tree_method)
    val (bucket_size)

    output:
    file("*.sp") 

    script:
    """
    ## Sum-of-Pairs Score ##
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.sp"
       
    """
}

process EASEL_INFO {
    tag "EASEL_INFO on $id"
    publishDir "${params.outdir}/easel"

    input:
    file (test_alignment)
    file(ref_alignment)
    val (align_method)
    val (tree_method)
    val (bucket_size)

    output:
    set file("${id}.${align_method}.${tree_method}.DPA.easel_INFO")
    file("${id}.${align_method}.${tree_method}.DPA.easel_AVG")
    file("${id}.${align_method}.${tree_method}.DPA.dump")

     shell:
     '''
     esl-alistat --icinfo !{id}.!{align_method}.!{tree_method}.DPA.easel_INFO !{id}.!{align_method}.!{tree_method}.DPA.dump

     awk 'NR > 8 && $1 !~/\\// { sum+= $3 } END {print "SUM: "sum"\\nAVG: "sum/(NR-9)}' !{id}.!{align_method}.!{tree_method}.DPA.easel_INFO > !{id}.!{align_method}.!{tree_method}.DPA.easel_AVG
     
     ## the first && is to skip first lines and the last one. The AVG is done -8 all the time execpt for the END print to "erase" the last "//" too.
     '''
}

process HOMOPLASY_TCOFFEE {
    tag "HOMOPLASY_TCOFFEE on $id"
    publishDir "${params.outdir}/homplasy"

    input:

    output:

     shell:
     '''
     '''
}