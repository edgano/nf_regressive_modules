#!/bin/bash nextflow
params.outdir = 'results'

include PRECOMPUTE_BLAST     from './preprocess.nf'   
include TCOFFEE_ALIGNER     from './generateAlignment.nf'   
include EVAL_ALIGNMENT      from './evaluateAlignment.nf'  
include EASEL_INFO          from './evaluateAlignment.nf'  
include GAPS_PROGRESSIVE    from './evaluateAlignment.nf'  
include METRICS             from './evaluateAlignment.nf' 

workflow PROG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
     
  main: 
    PROG_ALIGNER (seqs_and_trees, align_method)
   
    if (params.evaluate){
      refs_ch
        .cross (PROG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }
        
      EVAL_ALIGNMENT ("progressive", alignment_and_ref, PROG_ALIGNER.out.alignMethod, PROG_ALIGNER.out.treeMethod,"NA")
      EVAL_ALIGNMENT.out.tcScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.progressive.tcScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.spScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.progressive.spScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.colScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.progressive.colScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
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

workflow TCOFFEE_ANALYSIS {
  take:
    seqs
    tc_mode
     
  main: 
    PRECOMPUTE_BLAST (seqs)
    TCOFFEE_ALIGNER (seqs, tc_mode, PRECOMPUTE_BLAST.out.id)
}