#!/bin/bash nextflow
params.outdir = 'results'

include COMBINE_SEQS   from './preprocess.nf'    
include TREE_GENERATION   from './treeGeneration.nf'    
include REG_ALIGNER   from './generateAlignment.nf'   
include PROG_ALIGNER   from './generateAlignment.nf'   
include EVAL_ALIGNMENT   from './evaluateAlignment.nf'  

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
      TREE_GENERATION (seqs_ch, tree_methods) }
      else{/*define TREE_GENERATION.out*/}
    REG_ALIGNER (seqs_ch, align_methods, tree_methods, bucket_size, TREE_GENERATION.out)
    EVAL_ALIGNMENT (REG_ALIGNER.out.id, REG_ALIGNER.out.alignment, refs_ch, align_methods, tree_methods, bucket_size)
    //EASEL_INFO

  //emit: 
}