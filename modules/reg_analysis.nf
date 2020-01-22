#!/bin/bash nextflow
params.outdir = 'results'

include './combineSeqs' params(params)
include './treeGeneration' params(params)
include './generateAlignment' params(params)
include './evaluateAlignment' params(params)

workflow REG_ANALYSIS {
  get:
    seqs_ch
    refs_ch
    align_methods
    tree_methods
    bucket_size
     
  main: 
    //COMBINE_SEQS(seqs_ch, refs_ch) // need to combine seqs and ref by ID
    if (params.trees){
      TREE_GENERATION (seqs_ch, tree_methods)                                               }else{/*define TREE_GENERATION.out*/}
    if (params.regressive_align){
      REG_ALIGNER (seqs_ch, align_methods, tree_methods, bucket_size, TREE_GENERATION.out)  }else{/*define REG_ALIGNER.out*/}
      PROG_ALIGNER (seqs_ch, align_methods, tree_methods, TREE_GENERATION.out)
    //EVAL_ALIGNMENT (REG_ALIGNER.out, refs_ch, align_methods, tree_methods, bucket_size)

    //EASEL_INFO

  //emit: 
}