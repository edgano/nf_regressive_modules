#!/bin/bash nextflow

include TREE_GENERATION   from './treeGeneration.nf'    
include REG_ALIGNER       from './generateAlignment.nf'   
include PROG_ALIGNER      from './generateAlignment.nf'   
include EVAL_ALIGNMENT    from './evaluateAlignment.nf'  
include EASEL_INFO        from './evaluateAlignment.nf'  
include HOMOPLASY         from './evaluateAlignment.nf'  
include METRICS           from './evaluateAlignment.nf'  

workflow NAT_BIOTECH_ANALYSIS {
  take:
    seqs_ch
    refs_ch
    align_method
    tree_method
    bucket_size
    trees
     
  main: 
    if (!params.trees){
      TREE_GENERATION (seqs_ch, tree_method) 
    
      seqs_ch
          .cross(TREE_GENERATION.out)
          .map { it -> [ it[1][0], it[1][1], it[0][1], it[1][2] ] }
          .set { seqs_and_trees }
    }else{
      seqs_ch
        .cross(trees)
        .map { it -> [ it[1][0], it[1][1], it[0][1], it[1][2] ] }
        .set { seqs_and_trees }
    }
    REG_ALIGNER (seqs_and_trees, align_method, bucket_size)
   
    if (params.evaluate){
      refs_ch
        .cross (REG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      EVAL_ALIGNMENT ("regressive", alignment_and_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)
    }
    if (params.homoplasy){
      HOMOPLASY("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.metricFile)
    }
    if (params.easel){
      EASEL_INFO ("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)
    }
}
//REG_PROTOCOLS_ANALYSIS