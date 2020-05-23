#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'XXXXXX'.
 *
 *   XXXXXX is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   XXXXXX is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with XXXXXX.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main XXX pipeline script
 *
 * @authors
 * Edgar Garriga
 * Jose Espinosa-Carrasco
 */

//  example         https://github.com/nextflow-io/rnaseq-nf/tree/modules
/* 
 * enables modules 
 */
nextflow.preview.dsl = 2

/*
 * defaults parameter definitions
 */

//    ## subdatsets
//seq2improve="cryst,blmb,rrm,subt,ghf5,sdr,tRNA-synt_2b,zf-CCHH,egf,Acetyltransf,ghf13,p450,Rhodanese,aat,az,cytb,proteasome,GEL"
//params.seqs ="/users/cn/egarriga/datasets/homfam/combinedSeqs/{${seq2improve}}.fa"
top20fam=""

// input sequences to align in fasta format
//params.seqs = "$baseDir/data/*.fa"
params.seqs = 'https://raw.githubusercontent.com/edgano/datasets-test/homfam/seatoxin.fa'

//params.refs = "$baseDir/data/*.ref"
params.refs = 'https://raw.githubusercontent.com/edgano/datasets-test/homfam/seatoxin.ref'

//params.trees ="/Users/edgargarriga/CBCRG/nf_regressive_modules/results/trees/seatoxin.MBED.dnd"
params.trees = false
                      //CLUSTALO,FAMSA,MAFFT-FFTNS1
params.align_methods = "MAFFT-FFTNS1" 
                      //DPPARTTREE0,FAMSA-SLINK,MBED,PARTTREE
params.tree_methods = "FAMSA-SLINK"      

params.buckets = "10"

//  ## SLAVE parameters
                          //need to be lowercase -> direct to tcoffee
                          //mbed,parttree,famsadnd
params.slave_tree_methods="mbed,parttree,famsadnd" 

//  ## DYNAMIC parameters
params.dynamicX = "10000"
          //TODO -> make 2 list? one with aligners and the other with sizes?
          //      params.dynamicAlnList=["psicoffee_msa", "clustalo_msa", "famsa_msa"]
          //      params.dynamicMaxNseqList=[20, 10000, 1000000]
params.dynamicMasterAln="psicoffee_msa"
params.dynamicMasterSize="50"
params.dynamicSlaveAln="famsa_msa"
params.dynamicSlaveSize="100000000"
params.dynamicConfig=false

          //uniref50, pdb or path
params.db = "pdb"        

params.progressive_align = false
params.regressive_align = false 
params.pool_align=false        //<< TODO <- fix MAFFT on pool when nreg > datasetSize
params.slave_align=false
params.dynamic_align=true

params.evaluate=false
params.homoplasy=false
params.gapCount=false
params.metrics=false
params.easel=false

// output directory
params.outdir = "$baseDir/results"

// define database path
uniref_path = "/users/cn/egarriga/datasets/db/uniref50.fasta"   // cluster path
pdb_path = "/database/pdb/pdb_seqres.txt"                       // docker path

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
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (Aligned FASTA))              : ${params.refs}
         Input trees (NEWICK)                           : ${params.trees}
         Alignment methods                              : ${params.align_methods}
         Tree methods                                   : ${params.tree_methods}
         Bucket size                                    : ${params.buckets}
         --##--
         Generate Progressive alignments                : ${params.progressive_align}
         Generate Regressive alignments                 : ${params.regressive_align}
         Generate Slave tree alignments                 : ${params.slave_align}
                  Slave tree methods                    : ${params.slave_tree_methods}
         Generate Dynamic alignments                    : ${params.dynamic_align}
                  Dynamic size                          : ${params.dynamicX}
                  Dynamic config file                   : ${params.dynamicConfig}
                          master align - boundary       : ${params.dynamicMasterAln} - ${params.dynamicMasterSize}
                          slave align  - boundary       : ${params.dynamicSlaveAln} - ${params.dynamicSlaveSize}
                  Dynamic DDBB                          : ${params.db}
                  DDBB path                             : ${params.database_path}
         Generate Pool alignments                       : ${params.pool_align}
         --##--
         Perform evaluation? Requires reference         : ${params.evaluate}
         Check homoplasy? Only for regressive           : ${params.homoplasy}
         Check gapCount? For progressive                : ${params.gapCount}
         Check metrics?                                 : ${params.metrics}
         Check easel info?                              : ${params.easel}
         --##--
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()

// import analysis pipelines
include REG_ANALYSIS from './modules/reg_analysis'        params(params)
include PROG_ANALYSIS from './modules/prog_analysis'      params(params)
include SLAVE_ANALYSIS from './modules/reg_analysis'      params(params)
include DYNAMIC_ANALYSIS from './modules/reg_analysis'    params(params)
include POOL_ANALYSIS from './modules/reg_analysis'       params(params)

// Channels containing sequences
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

if ( params.refs ) {
  refs_ch = Channel.fromPath( params.refs ).map { item -> [ item.baseName, item] }
}

// Channels for user provided trees or empty channel if trees are to be generated [OPTIONAL]
if ( params.trees ) {
  trees = Channel.fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
}else { 
  Channel.empty().set { trees }
}

// tokenize params 
tree_method = params.tree_methods.tokenize(',')
align_method = params.align_methods.tokenize(',')
bucket_list = params.buckets.tokenize(',')
slave_method = params.slave_tree_methods.tokenize(',')
dynamicX = params.dynamicX   //TODO <- tokenize it ??

/* 
 * main script flow
 */
workflow pipeline {
    if (params.regressive_align){
      REG_ANALYSIS(seqs_ch, refs_ch, align_method, tree_method, bucket_list, trees)
    }
    if (params.progressive_align){
      PROG_ANALYSIS(seqs_ch, refs_ch, align_method, tree_method, trees)
    }
    if (params.slave_align){
      SLAVE_ANALYSIS(seqs_ch, refs_ch, align_method, tree_method, bucket_list, trees, slave_method)
    }
    if (params.dynamic_align){
      DYNAMIC_ANALYSIS(seqs_ch, refs_ch, tree_method, bucket_list, dynamicX, trees)
    }
    if (params.pool_align){
      POOL_ANALYSIS(seqs_ch, refs_ch, align_method, tree_method, bucket_list, trees)
    }
}

workflow {
  pipeline()
}

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
  //TODO script to generate CSV from individual files
}
