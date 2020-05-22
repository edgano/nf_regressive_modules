#!/bin/bash nextflow
params.outdir = 'results'

include COMBINE_SEQS        from './preprocess.nf'    
include TREE_GENERATION     from './treeGeneration.nf'    
include PROG_ALIGNER        from './generateAlignment.nf'   
include EVAL_ALIGNMENT      from './evaluateAlignment.nf'  
include EASEL_INFO          from './evaluateAlignment.nf'  
include GAPS_PROGRESSIVE    from './evaluateAlignment.nf'  
include METRICS             from './evaluateAlignment.nf' 

workflow PROG_ANALYSIS {
  take:
    seqs_ch
    refs_ch
    align_method
    tree_method
    trees
     
  main: 
    //COMBINE_SEQS(seqs_ch, refs_ch) // need to combine seqs and ref by ID
    if (!params.trees){
      TREE_GENERATION (seqs_ch, tree_method) 
        
      PROG_ALIGNER (seqs_ch, align_method, TREE_GENERATION.out.treeMethod, TREE_GENERATION.out.guideTree)
    }else{
      PROG_ALIGNER (seqs_ch, align_method, tree_method, trees)
    }

    refs_ch
        .cross (PROG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_plus_ref }

    if (params.evaluate){
      EVAL_ALIGNMENT ("progressive", alignment_plus_ref, PROG_ALIGNER.out.alignMethod, PROG_ALIGNER.out.treeMethod,"NA")
    }
    if (params.gapCount){
      GAPS_PROGRESSIVE("progressive", PROG_ALIGNER.out.alignmentFile, PROG_ALIGNER.out.alignMethod, PROG_ALIGNER.out.treeMethod,"NA")
    }
    if (params.metrics){
      METRICS("progressive", PROG_ALIGNER.out.alignmentFile, PROG_ALIGNER.out.alignMethod, PROG_ALIGNER.out.treeMethod,"NA", PROG_ALIGNER.out.metricFile)
    }
    if (params.easel){
      EASEL_INFO ("progressive", PROG_ALIGNER.out.alignmentFile, PROG_ALIGNER.out.alignMethod, PROG_ALIGNER.out.treeMethod,"NA")
    }
}