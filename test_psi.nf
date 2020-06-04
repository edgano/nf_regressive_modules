#!/usr/bin/env nextflow

/* 
 * enables modules 
 */
nextflow.preview.dsl = 2

/*
 * defaults parameter definitions
 */

// input sequences to align in fasta format
params.seqs = "${baseDir}/test/*.fa"

          //uniref50, pdb or path
params.db = "uniref50"        
          // define default database path
uniref_path = "/users/cn/egarriga/datasets/db/uniref50.fasta"   // cluster path
pdb_path    = "/database/pdb/pdb_seqres.txt"                       // docker path

params.numThreads = 1

// output directory
params.outdir = "$baseDir/results_test"

if (params.db=='uniref50'){
  params.database_path = uniref_path
}else if(params.db=='pdb'){
  params.database_path = pdb_path
}else{
  params.database_path = params.db
}

log.info """\
         PIPELINE  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)          : ${params.seqs}
         DDBB                             : ${params.db}
         DDBB path                        : ${params.database_path}
         numThreads                       : ${params.numThreads}
         --##--
         Output directory (DIRECTORY)     : ${params.outdir}
         """
         .stripIndent()

// import analysis pipelines
include BLASTP from './modules/bioinfoCommands'   params(params)

// Channels containing sequences
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

/*    main script flow  */
workflow pipeline {
  BLASTP(seqs_ch, params.numThreads)
}

workflow {
  pipeline()
}

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
