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
        trees_ch = TREE_GENERATION.out
    }else{
        trees_ch = trees
    }
    
    REG_ALIGNER (seqs_ch, align_method, bucket_size, trees_ch)

    if (params.evaluate){
      EVAL_ALIGNMENT ("regressive",REG_ALIGNER.out.id, REG_ALIGNER.out.alignmentFile, refs_ch, align_method, tree_method, bucket_size) 
      
      tcScore_csv = EVAL_ALIGNMENT.out.tcScore
                    .collectFile(name: "_tc.csv", newLine: true).view()
    }
    if (params.homoplasy){
      HOMOPLASY("regressive",REG_ALIGNER.out.id, REG_ALIGNER.out.alignmentFile, refs_ch, align_method, tree_method, bucket_size, REG_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("regressive",REG_ALIGNER.out.id, REG_ALIGNER.out.alignmentFile, refs_ch, align_method, tree_method, bucket_size, REG_ALIGNER.out.metricFile)
    }
    EASEL_INFO ("regressive",REG_ALIGNER.out.id, REG_ALIGNER.out.alignmentFile, refs_ch, align_method, tree_method, bucket_size)

    //emmit : tcScore_csv
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

    if (params.evaluate){
      EVAL_ALIGNMENT ("pool",POOL_ALIGNER.out.id, POOL_ALIGNER.out.alignmentFile, refs_ch, align_methods, tree_methods, bucket_size)
    }
    if (params.homoplasy){
      HOMOPLASY("pool",POOL_ALIGNER.out.id, POOL_ALIGNER.out.alignmentFile, refs_ch, align_methods, tree_methods, bucket_size, POOL_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("pool",POOL_ALIGNER.out.id, POOL_ALIGNER.out.alignmentFile, refs_ch, align_methods, tree_methods, bucket_size, POOL_ALIGNER.out.metricFile)
    }
    EASEL_INFO ("pool",POOL_ALIGNER.out.id, POOL_ALIGNER.out.alignmentFile, refs_ch, align_methods, tree_methods, bucket_size)

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

    if (params.evaluate){
      EVAL_ALIGNMENT ("slave",SLAVE_ALIGNER.out.id, SLAVE_ALIGNER.out.alignmentFile, refs_ch, align_methods, SLAVE_ALIGNER.out.tree_method, bucket_size)
    }
    if (params.homoplasy){
      HOMOPLASY("slave",SLAVE_ALIGNER.out.id, SLAVE_ALIGNER.out.alignmentFile, refs_ch, align_methods, SLAVE_ALIGNER.out.tree_method, bucket_size, SLAVE_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("slave",SLAVE_ALIGNER.out.id, SLAVE_ALIGNER.out.alignmentFile, refs_ch, align_methods, SLAVE_ALIGNER.out.tree_method, bucket_size, SLAVE_ALIGNER.out.metricFile)
    }
    EASEL_INFO ("slave",SLAVE_ALIGNER.out.id, SLAVE_ALIGNER.out.alignmentFile, refs_ch, align_methods, SLAVE_ALIGNER.out.tree_method, bucket_size)

  //emit: 
}