#!/bin/bash nextflow

include REG_ALIGNER       from './generateAlignment.nf'   
include PROG_ALIGNER      from './generateAlignment.nf'   
include EVAL_ALIGNMENT    from './evaluateAlignment.nf'  
include EASEL_INFO        from './evaluateAlignment.nf'  
include HOMOPLASY         from './evaluateAlignment.nf'  
include METRICS           from './evaluateAlignment.nf'  

workflow NAT_BIOTECH_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
     
  main: 
    REG_ALIGNER (seqs_and_trees, align_method, bucket_size)
   
      refs_ch
        .cross (REG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      EVAL_ALIGNMENT ("regressive", alignment_and_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)

      HOMOPLASY("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.homoplasyFile)

      METRICS("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.metricFile)

      EASEL_INFO ("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)
    
}
workflow REG_PROTOCOLS_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
     
  main: 
}

workflow TCOFFEE_PROTOCOLS_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
     
  main: 
}