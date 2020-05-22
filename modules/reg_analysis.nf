#!/bin/bash nextflow
params.outdir = 'results_REG'

include COMBINE_SEQS      from './preprocess.nf'    
include TREE_GENERATION   from './treeGeneration.nf'    
include REG_ALIGNER       from './generateAlignment.nf'   
include EVAL_ALIGNMENT    from './evaluateAlignment.nf'  
include EASEL_INFO        from './evaluateAlignment.nf'  
include HOMOPLASY         from './evaluateAlignment.nf'  
include METRICS           from './evaluateAlignment.nf'  

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

    //REG_ALIGNER.out.alignmentFile.view()

    if (params.evaluate){
      EVAL_ALIGNMENT ("regressive", alignment_plus_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)
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

include SLAVE_ALIGNER   from './generateAlignment.nf'   
workflow SLAVE_ANALYSIS {
  take:
    seqs_ch
    refs_ch
    align_method
    tree_method
    bucket_size
    trees
    slave_method
     
  main: 

      //COMBINE_SEQS(seqs_ch, refs_ch) // need to combine seqs and ref by ID
    if (!params.trees){
      TREE_GENERATION (seqs_ch, tree_method) 
        
      SLAVE_ALIGNER (seqs_ch, align_method, bucket_size, TREE_GENERATION.out.treeMethod, TREE_GENERATION.out.guideTree, slave_method)
    }else{
      SLAVE_ALIGNER (seqs_ch, align_method, bucket_size, tree_method, trees, slave_method)
    }

    refs_ch
        .cross (SLAVE_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_plus_ref }

    if (params.evaluate){
      EVAL_ALIGNMENT ("slave", alignment_plus_ref, SLAVE_ALIGNER.out.alignMethod, SLAVE_ALIGNER.out.treeMethod, SLAVE_ALIGNER.out.bucketSize)
    }
    if (params.homoplasy){
      HOMOPLASY("slave", SLAVE_ALIGNER.out.alignmentFile, SLAVE_ALIGNER.out.alignMethod, SLAVE_ALIGNER.out.treeMethod, SLAVE_ALIGNER.out.bucketSize, SLAVE_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("slave", SLAVE_ALIGNER.out.alignmentFile, SLAVE_ALIGNER.out.alignMethod, SLAVE_ALIGNER.out.treeMethod, SLAVE_ALIGNER.out.bucketSize, SLAVE_ALIGNER.out.metricFile)
    }
    if (params.easel){
      EASEL_INFO ("slave", SLAVE_ALIGNER.out.alignmentFile, SLAVE_ALIGNER.out.alignMethod, SLAVE_ALIGNER.out.treeMethod, SLAVE_ALIGNER.out.bucketSize)
    }
}

include GENERATE_DYNAMIC_CONFIG      from './preprocess.nf'    
include DYNAMIC_ALIGNER             from './generateAlignment.nf'   
workflow DYNAMIC_ANALYSIS {
  take:
    seqs_ch
    refs_ch
    tree_method
    bucket_size
    dynamicX
    trees
     
  main: 
      //COMBINE_SEQS(seqs_ch, refs_ch) // need to combine seqs and ref by ID

    if(params.dynamicConfig){
      GENERATE_DYNAMIC_CONFIG(params.dynamicMasterAln, params.dynamicMasterSize, params.dynamicSlaveAln, params.dynamicSlaveSize)
      align_method="CONFIG"
      configFile = GENERATE_DYNAMIC_CONFIG.out.configFile
      configValues = GENERATE_DYNAMIC_CONFIG.out.configValues
      dynamicValues = "${params.dynamicMasterAln}.${params.dynamicMasterSize}_${params.dynamicSlaveAln}.${params.dynamicSlaveSize}"
    }else{
      align_method="DEFAULT"
      configFile = "/"
      configValues=["","","",""]
      dynamicValues = "DEFAULT"
    }

    if (!params.trees){
      TREE_GENERATION (seqs_ch, tree_method) 
        
      DYNAMIC_ALIGNER (seqs_ch, align_method, bucket_size, dynamicX, TREE_GENERATION.out.treeMethod, TREE_GENERATION.out.guideTree, configFile, configValues,dynamicValues)
    }else{
      DYNAMIC_ALIGNER (seqs_ch, align_method, bucket_size, dynamicX, tree_method, trees, configFile, configValues,dynamicValues)
    }

    refs_ch
        .cross (DYNAMIC_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_plus_ref }

    if (params.evaluate){
      EVAL_ALIGNMENT ("dynamic", alignment_plus_ref, DYNAMIC_ALIGNER.out.alignMethod, DYNAMIC_ALIGNER.out.treeMethod, DYNAMIC_ALIGNER.out.bucketSize)
    }
    if (params.homoplasy){
      HOMOPLASY("dynamic", DYNAMIC_ALIGNER.out.alignmentFile, DYNAMIC_ALIGNER.out.alignMethod, DYNAMIC_ALIGNER.out.treeMethod, DYNAMIC_ALIGNER.out.bucketSize, DYNAMIC_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("dynamic", DYNAMIC_ALIGNER.out.alignmentFile, DYNAMIC_ALIGNER.out.alignMethod, DYNAMIC_ALIGNER.out.treeMethod, DYNAMIC_ALIGNER.out.bucketSize, DYNAMIC_ALIGNER.out.metricFile)
    }
    if (params.easel){
      EASEL_INFO ("dynamic", DYNAMIC_ALIGNER.out.alignmentFile, DYNAMIC_ALIGNER.out.alignMethod, DYNAMIC_ALIGNER.out.treeMethod, DYNAMIC_ALIGNER.out.bucketSize)
    }
}

include POOL_ALIGNER   from './generateAlignment.nf'   
workflow POOL_ANALYSIS {
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
        
      POOL_ALIGNER (seqs_ch, align_method, bucket_size, TREE_GENERATION.out.treeMethod, TREE_GENERATION.out.guideTree)
    }else{
      POOL_ALIGNER (seqs_ch, align_method, bucket_size, tree_method, trees)
    }

    refs_ch
        .cross (POOL_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_plus_ref }

    if (params.evaluate){
      EVAL_ALIGNMENT ("pool", alignment_plus_ref, POOL_ALIGNER.out.alignMethod, POOL_ALIGNER.out.treeMethod, POOL_ALIGNER.out.bucketSize)
    }
    if (params.homoplasy){
      HOMOPLASY("pool", POOL_ALIGNER.out.alignmentFile, POOL_ALIGNER.out.alignMethod, POOL_ALIGNER.out.treeMethod, POOL_ALIGNER.out.bucketSize, POOL_ALIGNER.out.homoplasyFile)
    }
    if (params.metrics){
      METRICS("pool", POOL_ALIGNER.out.alignmentFile, POOL_ALIGNER.out.alignMethod, POOL_ALIGNER.out.treeMethod, POOL_ALIGNER.out.bucketSize, POOL_ALIGNER.out.metricFile)
    }
    if (params.easel){
      EASEL_INFO ("pool", POOL_ALIGNER.out.alignmentFile, POOL_ALIGNER.out.alignMethod, POOL_ALIGNER.out.treeMethod, POOL_ALIGNER.out.bucketSize)
    }
}