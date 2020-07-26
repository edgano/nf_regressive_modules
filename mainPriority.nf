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
seq2improve="cryst,blmb,rrm,subt,ghf5,sdr,tRNA-synt_2b,zf-CCHH,egf,Acetyltransf,ghf13,p450,Rhodanese,aat,az,cytb,proteasome,GEL"
top20fam="gluts,myb_DNA-binding,tRNA-synt_2b,biotin_lipoyl,hom,ghf13,aldosered,hla,Rhodanese,PDZ,blmb,rhv,p450,adh,aat,rrm,Acetyltransf,sdr,zf-CCHH,rvp"
smallTest="seatoxin,hip,scorptoxin,cyt3,rnasemam,bowman,toxin,ghf11,TNF,sti"


params.seqs ="/users/cn/egarriga/datasets/homfam/combinedSeqs/{${smallTest}}.fa"

// input sequences to align in fasta format
//params.seqs = "/users/cn/egarriga/datasets/homfam/combinedSeqs/*.fa"

params.refs = "/users/cn/egarriga/datasets/homfam/refs/*.ref"

params.trees ="/users/cn/egarriga/datasets/homfam/trees/*.CLUSTALO.dnd"
//params.trees = false
                      //CLUSTALO,FAMSA,MAFFT-FFTNS1
params.align_methods = "PSI,PSI_FLAG"//,FAMSA,MAFFT-FFTNS1" 
                      //MAFFT-DPPARTTREE0,FAMSA-SLINK,MBED,MAFFT-PARTTREE
params.tree_methods = "MBED"      //TODO -> reuse trees for multiple methods.

params.buckets = "50"

//  ## SLAVE parameters
                          //need to be lowercase -> direct to tcoffee
                          //mbed,parttree,famsadnd
params.slave_tree_methods="mbed,famsadnd,parttree" 

//  ## DYNAMIC parameters
params.dynamicX = "100000"
          //TODO -> make 2 list? one with aligners and the other with sizes? (to have more than 2 aligners)
params.dynamicMasterAln="psicoffee_msa"
params.dynamicMasterSize="50"
params.dynamicSlaveAln="famsa_msa"
params.dynamicSlaveSize="100000000"
params.dynamicConfig=false

          //uniref50, pdb or path
params.db = "uniref50"        
params.blastOutdir = "$baseDir/Blast"

params.progressive_align = false
params.regressive_align = false          
params.pool_align=false                  
params.slave_align=false    
params.dynamic_align=true                

params.evaluate=true
params.homoplasy=true
params.gapCount=false
params.metrics=true
params.easel=false

// output directory
params.outdir = "$baseDir/resultsPRIORITY"

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
include TREE_GENERATION from './modules/treeGeneration'        params(params)
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
bucket_list = params.buckets.toString().tokenize(',') //int to string
slave_method = params.slave_tree_methods.tokenize(',')
dynamicX = params.dynamicX.toString().tokenize(',') //int to string

/* 
 * main script flow
 */
workflow pipeline {

    if (!params.trees){
      TREE_GENERATION (seqs_ch, tree_method) 
      seqs_ch
        .cross(TREE_GENERATION.out)
        .map { it -> [ it[1][0], it[1][1], it[0][1], it[1][2] ] }
        .set { seqs_and_trees }
    }else{
      seqs_ch
        .cross(trees)
        .map { it -> [ it[1][0], it[1][1], it[0][1], it[1][2] ] }
        .set { seqs_and_trees }
    }
    // REGRESSIVE psicoffee & psicoffee_flag
        REG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)

    // PROGRESSIVE tcoffee
        PROG_ANALYSIS(seqs_and_trees, refs_ch, "TCOFFEE", tree_method)

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
