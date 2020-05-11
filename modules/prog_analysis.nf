#!/bin/bash nextflow
params.outdir = 'results'

include COMBINE_SEQS   from './preprocess.nf'    
include TREE_GENERATION   from './treeGeneration.nf'    
include PROG_ALIGNER   from './generateAlignment.nf'   
include EVAL_ALIGNMENT   from './evaluateAlignment.nf'  
include EASEL_INFO   from './evaluateAlignment.nf'  

workflow PROG_ANALYSIS {
  take:
    seqs_ch
    refs_ch
    align_methods
    tree_methods
     
  main: 
    //COMBINE_SEQS(seqs_ch, refs_ch) // need to combine seqs and ref by ID
    if (params.trees){
      TREE_GENERATION (seqs_ch, tree_methods) }
      else{/*define TREE_GENERATION.out*/}
    
    PROG_ALIGNER (seqs_ch, align_methods, tree_methods, TREE_GENERATION.out)

    if (params.evaluate){
      EVAL_ALIGNMENT ("progressive",PROG_ALIGNER.out.id, PROG_ALIGNER.out.alignment, refs_ch, align_methods, tree_methods, "NA")
    }
    EASEL_INFO ("progressive",PROG_ALIGNER.out.id, PROG_ALIGNER.out.alignment, refs_ch, align_methods, tree_methods, "NA")

  //emit: 
}