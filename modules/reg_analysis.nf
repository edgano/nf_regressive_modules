#!/bin/bash nextflow
params.outdir = 'results'

include COMBINE_SEQS   from './preprocess.nf'    params(params)
include TREE_GENERATION   from './treeGeneration.nf'    params(params)
include REG_ALIGNER   from './generateAlignment.nf'    params(params)
include PROG_ALIGNER   from './generateAlignment.nf'    params(params)
include EVAL_ALIGNMENT   from './evaluateAlignment.nf'    params(params)

workflow REG_ANALYSIS {
  take:
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