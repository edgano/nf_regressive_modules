#!/bin/bash nextflow
params.outdir = 'results'

include COMBINE_SEQS   from './preprocess.nf'    
include TREE_GENERATION   from './treeGeneration.nf'    
include PROG_ALIGNER   from './generateAlignment.nf'   
include EVAL_ALIGNMENT   from './evaluateAlignment.nf'  
include GAPS_PROGRESSIVE   from './evaluateAlignment.nf'  
include METRICS   from './evaluateAlignment.nf'  
include EASEL_INFO   from './evaluateAlignment.nf'  

workflow PROG_ANALYSIS {
  take:
    seqs_ch
    refs_ch
    align_methods
    tree_methods
    trees
     
  main: 
    //COMBINE_SEQS(seqs_ch, refs_ch) // need to combine seqs and ref by ID
    if (!params.trees){
        TREE_GENERATION (seqs_ch, tree_methods) 
        trees_ch = TREE_GENERATION.out
    }else{
        trees_ch = trees
    }
    
    PROG_ALIGNER (seqs_ch, align_methods, trees_ch)

    refs_ch
        .cross (PROG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { prog_alignment_plus_ref }

    if (params.evaluate){
      EVAL_ALIGNMENT ("progressive", prog_alignment_plus_ref, align_methods, tree_methods, "NA")
    }
    if (params.gapCount){
      GAPS_PROGRESSIVE("progressive", PROG_ALIGNER.out.id, PROG_ALIGNER.out.alignmentFile, align_methods, tree_methods, "NA")
    }
    if (params.metrics){
      METRICS("progressive", prog_alignment_plus_ref, align_methods, tree_methods, "NA", PROG_ALIGNER.out.metricFile)
    }
    EASEL_INFO ("progressive", prog_alignment_plus_ref, align_methods, tree_methods, "NA")

  //emit: 
}