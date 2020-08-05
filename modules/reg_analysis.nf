#!/bin/bash nextflow
//params.outdir = 'results_REG'

include {COMBINE_SEQS}      from './preprocess.nf'    
include {REG_ALIGNER}       from './generateAlignment.nf'   
include {EVAL_ALIGNMENT}    from './modules_evaluateAlignment.nf'
include {EASEL_INFO}        from './modules_evaluateAlignment.nf'
include {HOMOPLASY}         from './modules_evaluateAlignment.nf'
include {METRICS}           from './modules_evaluateAlignment.nf'

workflow REG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
    
     
  main: 
    REG_ALIGNER (seqs_and_trees, align_method, bucket_size)
   
    if (params.evaluate){
      refs_ch
        .cross (REG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      EVAL_ALIGNMENT ("regressive", alignment_and_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)
      EVAL_ALIGNMENT.out.tcScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.tcScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.spScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.spScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.colScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.colScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
                    
    }
    if (params.homoplasy){
      HOMOPLASY("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.homoplasyFile)
      HOMOPLASY.out.homoFiles
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text};${it[10].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.homo.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")  
    }

    def metrics_regressive = params.metrics? METRICS("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.metricFile) : Channel.empty()
    if (params.metrics) {
        metrics_regressive.metricFiles
                          .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text}" }
                          .collectFile(name: "${workflow.runName}.regressive.metrics.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }

    def easel_info = params.easel? EASEL_INFO ("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize) : Channel.empty()
    if (params.easel) {
        easel_info.easelFiles
                  .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[6].text};${it[7].text}" }
                  .collectFile(name: "${workflow.runName}.regressive.easel.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }

    emit:
    alignment = REG_ALIGNER.out.alignmentFile
    metrics = metrics_regressive
    easel = easel_info

}

include {SLAVE_ALIGNER}   from './generateAlignment.nf'   
workflow SLAVE_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
    slave_method
     
  main: 
    SLAVE_ALIGNER (seqs_and_trees, align_method, bucket_size, slave_method)
   
    if (params.evaluate){
      refs_ch
        .cross (SLAVE_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      EVAL_ALIGNMENT ("slave", alignment_and_ref, SLAVE_ALIGNER.out.alignMethod, SLAVE_ALIGNER.out.treeMethod, SLAVE_ALIGNER.out.bucketSize)
      EVAL_ALIGNMENT.out.tcScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.slave.tcScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.spScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.slave.spScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.colScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.slave.colScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }
    if (params.homoplasy){
      HOMOPLASY("slave", SLAVE_ALIGNER.out.alignmentFile, SLAVE_ALIGNER.out.alignMethod, SLAVE_ALIGNER.out.treeMethod, SLAVE_ALIGNER.out.bucketSize, SLAVE_ALIGNER.out.homoplasyFile)
      HOMOPLASY.out.homoFiles
                  .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text};${it[10].text}" }
                  .collectFile(name: "${workflow.runName}.slave.homo.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")     
    }
    if (params.metrics){
      METRICS("slave", SLAVE_ALIGNER.out.alignmentFile, SLAVE_ALIGNER.out.alignMethod, SLAVE_ALIGNER.out.treeMethod, SLAVE_ALIGNER.out.bucketSize, SLAVE_ALIGNER.out.metricFile)
      METRICS.out.metricFiles
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text}" }
                    .collectFile(name: "${workflow.runName}.slave.metrics.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }
    if (params.easel){
      EASEL_INFO ("slave", SLAVE_ALIGNER.out.alignmentFile, SLAVE_ALIGNER.out.alignMethod, SLAVE_ALIGNER.out.treeMethod, SLAVE_ALIGNER.out.bucketSize)
      EASEL_INFO.out.easelFiles
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[6].text};${it[7].text}" }
                    .collectFile(name: "${workflow.runName}.slave.easel.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }

  emit:
  alignment = SLAVE_ALIGNER.out.alignmentFile
}

include {GENERATE_DYNAMIC_CONFIG}      from './preprocess.nf'    
include {DYNAMIC_ALIGNER}             from './generateAlignment.nf'   
workflow DYNAMIC_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    tree_method
    bucket_size
    dynamicX
     
  main: 
// ######
//      TODO <- refactor this mess
// ######
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

    DYNAMIC_ALIGNER (seqs_and_trees, align_method, bucket_size, dynamicX, configFile, configValues, dynamicValues)
    
    if (params.evaluate){
      refs_ch
        .cross (DYNAMIC_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      EVAL_ALIGNMENT ("dynamic", alignment_and_ref, DYNAMIC_ALIGNER.out.alignMethod, DYNAMIC_ALIGNER.out.treeMethod, DYNAMIC_ALIGNER.out.bucketSize)
      EVAL_ALIGNMENT.out.tcScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.dynamic.tcScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.spScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.dynamic.spScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.colScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.dynamic.colScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }
    if (params.homoplasy){
      HOMOPLASY("dynamic", DYNAMIC_ALIGNER.out.alignmentFile, DYNAMIC_ALIGNER.out.alignMethod, DYNAMIC_ALIGNER.out.treeMethod, DYNAMIC_ALIGNER.out.bucketSize, DYNAMIC_ALIGNER.out.homoplasyFile)
      HOMOPLASY.out.homoFiles
                  .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text};${it[10].text}" }
                  .collectFile(name: "${workflow.runName}.dynamic.homo.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")  
    }
    if (params.metrics){
      METRICS("dynamic", DYNAMIC_ALIGNER.out.alignmentFile, DYNAMIC_ALIGNER.out.alignMethod, DYNAMIC_ALIGNER.out.treeMethod, DYNAMIC_ALIGNER.out.bucketSize, DYNAMIC_ALIGNER.out.metricFile)
      METRICS.out.metricFiles
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text}" }
                    .collectFile(name: "${workflow.runName}.dynamic.metrics.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }
    if (params.easel){
      EASEL_INFO ("dynamic", DYNAMIC_ALIGNER.out.alignmentFile, DYNAMIC_ALIGNER.out.alignMethod, DYNAMIC_ALIGNER.out.treeMethod, DYNAMIC_ALIGNER.out.bucketSize)
      EASEL_INFO.out.easelFiles
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[6].text};${it[7].text}" }
                    .collectFile(name: "${workflow.runName}.dynamic.easel.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")    
    }

  emit:
  alignment = DYNAMIC_ALIGNER.out.alignmentFile
}

include {POOL_ALIGNER}   from './generateAlignment.nf'   
workflow POOL_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
     
  main: 
    POOL_ALIGNER (seqs_and_trees, align_method, bucket_size)
   
    if (params.evaluate){
      refs_ch
        .cross (POOL_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }
    
      EVAL_ALIGNMENT ("pool", alignment_and_ref, POOL_ALIGNER.out.alignMethod, POOL_ALIGNER.out.treeMethod, POOL_ALIGNER.out.bucketSize)
      EVAL_ALIGNMENT.out.tcScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.pool.tcScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.spScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.pool.spScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.colScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.pool.colScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }
    if (params.homoplasy){
      HOMOPLASY("pool", POOL_ALIGNER.out.alignmentFile, POOL_ALIGNER.out.alignMethod, POOL_ALIGNER.out.treeMethod, POOL_ALIGNER.out.bucketSize, POOL_ALIGNER.out.homoplasyFile)
      HOMOPLASY.out.homoFiles
                  .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text};${it[10].text}" }
                  .collectFile(name: "${workflow.runName}.pool.homo.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")  
    }
    if (params.metrics){
      METRICS("pool", POOL_ALIGNER.out.alignmentFile, POOL_ALIGNER.out.alignMethod, POOL_ALIGNER.out.treeMethod, POOL_ALIGNER.out.bucketSize, POOL_ALIGNER.out.metricFile)
      METRICS.out.metricFiles
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text}" }
                    .collectFile(name: "${workflow.runName}.pool.metrics.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }
    if (params.easel){
      EASEL_INFO ("pool", POOL_ALIGNER.out.alignmentFile, POOL_ALIGNER.out.alignMethod, POOL_ALIGNER.out.treeMethod, POOL_ALIGNER.out.bucketSize)
      EASEL_INFO.out.easelFiles
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[6].text};${it[7].text}" }
                    .collectFile(name: "${workflow.runName}.pool.easel.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")    
    }

  emit:
  alignment = POOL_ALIGNER.out.alignmentFile
}