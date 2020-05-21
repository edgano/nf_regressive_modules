#!/bin/bash nextflow
params.outdir = 'results_REG'

include COMBINE_SEQS   from './preprocess.nf'    
include TREE_GENERATION   from './treeGeneration.nf'    
include REG_ALIGNER   from './generateAlignment.nf'   
include EVAL_ALIGNMENT   from './evaluateAlignment.nf'  
include EASEL_INFO   from './evaluateAlignment.nf'  
include HOMOPLASY   from './evaluateAlignment.nf'  
include METRICS   from './evaluateAlignment.nf'  

workflow REG_ANALYSIS {
  take:
    seqs_ch
    refs_ch
    align_method
    tree_method
    bucket_size
    trees
     
  main: 
    //COMBINE_SEQS(seqs_ch, refs_ch) // need to combine seqs and ref by ID
    if (!params.trees){
      TREE_GENERATION (seqs_ch, tree_method) 
        
      REG_ALIGNER (seqs_ch, align_method, bucket_size, TREE_GENERATION.out.treeMethod, TREE_GENERATION.out.guideTree)
    }else{

      REG_ALIGNER (seqs_ch, align_method, bucket_size, tree_method, trees)
    }

    refs_ch
        .cross (REG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_plus_ref }

    //REG_ALIGNER.out.alignmentFile.cross(refs_ch).view()

    if (params.evaluate){
      EVAL_ALIGNMENT ("regressive", alignment_plus_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)
    }
    return
    if (params.homoplasy){
      HOMOPLASY("regressive", alignment_plus_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("regressive", alignment_plus_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.metricFile)
    }
    // EASEL_INFO ("regressive", alignment_plus_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)

}

include POOL_ALIGNER   from './generateAlignment.nf'   
workflow POOL_ANALYSIS {
  take:
    seqs_ch
    refs_ch
    align_methods
    tree_methods
    bucket_size
    trees
     
  main: 
    //COMBINE_SEQS(seqs_ch, refs_ch) // need to combine seqs and ref by ID
    if (!params.trees){
        TREE_GENERATION (seqs_ch, tree_methods) 
        trees_ch = TREE_GENERATION.out
    }else{
        trees_ch = trees
    }
    
    POOL_ALIGNER (seqs_ch, align_methods, bucket_size, trees_ch)

    refs_ch
        .cross (POOL_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { pool_alignment_plus_ref }

    if (params.evaluate){
      EVAL_ALIGNMENT ("pool", pool_alignment_plus_ref, align_methods, tree_methods, bucket_size)
    }
    if (params.homoplasy){
      HOMOPLASY ("pool", pool_alignment_plus_ref, align_methods, tree_methods, bucket_size, POOL_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("pool", pool_alignment_plus_ref, align_methods, tree_methods, bucket_size, POOL_ALIGNER.out.metricFile)
    }
    EASEL_INFO ("pool", pool_alignment_plus_ref, align_methods, tree_methods, bucket_size)

  //emit: 
}

include SLAVE_ALIGNER   from './generateAlignment.nf'   
workflow SLAVE_ANALYSIS {
  take:
    seqs_ch
    refs_ch
    align_methods
    tree_methods
    bucket_size
    trees
    slave_method
     
  main: 
    //COMBINE_SEQS(seqs_ch, refs_ch) // need to combine seqs and ref by ID
    if (!params.trees){
        TREE_GENERATION (seqs_ch, tree_methods) 
        trees_ch = TREE_GENERATION.out
    }else{
        trees_ch = trees
    }
    
    SLAVE_ALIGNER (seqs_ch, align_methods, bucket_size, trees_ch, slave_method)

    refs_ch
        .cross (SLAVE_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { slave_alignment_plus_ref }

    if (params.evaluate){
      EVAL_ALIGNMENT ("slave", slave_alignment_plus_ref, align_methods, SLAVE_ALIGNER.out.tree_method, bucket_size)
    }
    if (params.homoplasy){
      HOMOPLASY("slave", slave_alignment_plus_ref, align_methods, SLAVE_ALIGNER.out.tree_method, bucket_size, SLAVE_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("slave", slave_alignment_plus_ref, align_methods, SLAVE_ALIGNER.out.tree_method, bucket_size, SLAVE_ALIGNER.out.metricFile)
    }
    EASEL_INFO ("slave", slave_alignment_plus_ref, align_methods, SLAVE_ALIGNER.out.tree_method, bucket_size)

  //emit: 
}