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
//params.seqs ="/users/cn/egarriga/datasets/homfam/combinedSeqs/{${seq2improve}}.fa"

// input sequences to align in fasta format
params.seqs = "/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/sh3.fasta"

params.refs = "/users/cn/egarriga/datasets/homfam/refs/*.ref"

params.trees ="/users/cn/egarriga/datasets/homfam/trees/*.{FAMSA,CLUSTALO,MAFFT_PARTTREE}.dnd"
//params.trees = false
                      //TODO FIX -> reg_UPP
                      //CLUSTALO,FAMSA,MAFFT-FFTNS1,MAFFT-GINSI,MAFFT-SPARSECORE,MAFFT,MSAPROBS,PROBCONS,TCOFFEE,UPP,MUSCLE
<<<<<<< HEAD
params.align_methods = "CLUSTALO"
                      
=======
params.align_methods = "CLUSTALO,FAMSA,MAFFT-FFTNS1"   

>>>>>>> bebf4fa34f11a60205ccf1ebaf050708dfebd9de
//CLUSTALW-QUICK,CLUSTALW                    
//FAMSA-SLINK,FAMSA-SLINKmedoid,FAMSA-SLINKparttree,FAMSA-UPGMA,FAMSA-UPGMAmedoid,FAMSA-UPGMAparttree   
//MAFFT-DPPARTTREE0,MAFFT-DPPARTTREE1,MAFFT-DPPARTTREE2,MAFFT-DPPARTTREE2size
//MAFFT-FASTAPARTTREE,MAFFT-FFTNS1,MAFFT-FFTNS1mem,MAFFT-FFTNS2,MAFFT-FFTNS2mem
//MAFFT-PARTTREE0,MAFFT-PARTTREE1,MAFFT-PARTTREE2,MAFFT-PARTTREE2size
//MAFFT,MBED
//TCOFFEE-BLENGTH,TCOFFEE-ISWLCAT,TCOFFEE-KM,TCOFFEE-LONGCAT,TCOFFEE-NJ,TCOFFEE-REG,TCOFFEE-SHORTCAT,TCOFFEE-SWL,TCOFFEE-SWLcat,TCOFFEE-UPGMA

//TODO -> test tcoffee trees
//                    CLUSTALW-QUICK,CLUSTALW  -> not working on PROG bc they are not rooted
                      //MAFFT-DPPARTTREE0,FAMSA-SLINK,MBED,MAFFT-PARTTREE0
<<<<<<< HEAD
params.tree_methods = "MBED"      
=======
params.tree_methods = "FAMSA-SLINK"      
>>>>>>> bebf4fa34f11a60205ccf1ebaf050708dfebd9de

params.buckets = "1000"

//  ## SLAVE parameters               //need to be lowercase -> direct to tcoffee
                          //mbed,parttree,famsadnd,cwdnd,dpparttree,fastparttree,mafftdnd,fftns1dnd,fftns2dnd,nj 
                          //mbed,parttree,famsadnd
params.slave_tree_methods="mbed,parttree,famsadnd" 

//  ## DYNAMIC parameters
params.dynamicX = "10000"                 //TODO -> make 2 list? one with aligners and the other with sizes? (to have more than 2 aligners)
params.dynamicMasterAln="psicoffee_msa"
params.dynamicMasterSize="50"
params.dynamicSlaveAln="famsa_msa"
params.dynamicSlaveSize="100000000"
params.dynamicConfig=true


// ## TCOFFEE
                  //3DALIGN,3DCOFFEE,3DMALIGN,ACCURATE,DEFAULT,EXPRESSO,FMCOFFEE,MCOFFEE,PROCOFFEE,PSICOFFEE,QUICKALN,RCOFFEE_CONSAN,RCOFFEE,TRMSD"
  //need template -> 3DMALIGN
params.tc_modes = "DEFAULT"
params.templates = '/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/sh3.template_file'
params.pdb = '/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/PDBs/*.pdb'
params.libs = ''
params.pairFile = ''
params.params4tcoffee = ''   
params.cache_path = ''

          //uniref50, pdb or path
params.db = "pdb"        
          // define default database path
uniref_path = "/users/cn/egarriga/datasets/db/uniref50.fasta"   // cluster path
pdb_path = "/database/pdb/pdb_seqres.txt"                       // docker path


params.progressive_align = false
<<<<<<< HEAD
params.regressive_align = true           
params.pool_align=false                  
params.slave_align=false   
params.dynamic_align=false               
=======
params.regressive_align = false           
params.pool_align=false                  
params.slave_align=false   
params.dynamic_align=false   
params.tcoffee_align = true            
>>>>>>> bebf4fa34f11a60205ccf1ebaf050708dfebd9de

params.evaluate=true
params.homoplasy=true
params.gapCount=true
params.metrics=true
params.easel=true

// output directory
params.outdir = "$baseDir/results"

//blast call cached
params.generateBlast=false
params.blastOutdir="$baseDir/blast"

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
         Generate TCoffee alignments                    : ${params.tcoffee_align}
                  TCoffee modes                         : ${params.tc_modes}
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
//include REG_ANALYSIS from './modules/reg_analysis'        params(params)
include REG_ANALYSIS from './modules/regressiveAnalysis'        params(params)
include PROG_ANALYSIS from './modules/prog_analysis'      params(params)
include SLAVE_ANALYSIS from './modules/reg_analysis'      params(params)
include DYNAMIC_ANALYSIS from './modules/reg_analysis'    params(params)
include POOL_ANALYSIS from './modules/reg_analysis'       params(params)
include TCOFFEE_CI from './modules/prog_analysis'      params(params)

seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

if ( params.refs ) {
  refs_ch = Channel.fromPath( params.refs ).map { item -> [ item.baseName, item] }
}

if ( params.templates ) {
  Channel
  .fromPath(params.templates, checkIfExists: true)
  .map { item -> [ item.simpleName, item] }
  .set { templates_ch}
}else { 
  Channel.empty()
  .set { templates_ch }
}
if ( params.pdb ) {
  Channel
  .fromPath(params.pdb)
  .map { item -> [ item.simpleName , item] }
  .groupTuple()
  .set { pdbFiles}
}else { 
  Channel.empty()
    .set { pdbFiles }
}
if ( params.pairFile ) {
  Channel
  .fromPath(params.pairFile)
  .map { item -> [ item.simpleName , item] }
  .groupTuple()
  .set { pair_file}
}else { 
  Channel.empty()
    .set { pair_file }
}
if ( params.libs ) {
  Channel
  .fromPath(params.libs)
  .map { item -> [ item.simpleName , item] }
  .groupTuple()
  .set { libs_ch }
}else { 
  Channel.empty()
    .set { libs_ch }
}

if ( params.trees ) {
  trees_ch = Channel.fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
}else { 
  Channel.empty().set { trees_ch }
}

// tokenize params 
tree_method = params.tree_methods.tokenize(',')
align_method = params.align_methods.tokenize(',')
bucket_list = params.buckets.toString().tokenize(',')     //int to string
slave_method = params.slave_tree_methods.tokenize(',')
dynamicX = params.dynamicX.toString().tokenize(',')       //int to string

tc_mode = params.tc_modes.tokenize(',')

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
    if (params.regressive_align){
      REG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)
    }
    if (params.progressive_align){
      PROG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method)
    }
    if (params.slave_align){
      SLAVE_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list, slave_method)
    }
    if (params.dynamic_align){
      DYNAMIC_ANALYSIS(seqs_and_trees, refs_ch, tree_method, bucket_list, dynamicX)
    }
    if (params.pool_align){
      POOL_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)
    }
    if (params.tcoffee_align){

      seqs_and_trees
        .join(templates_ch, remainder: true)
        .set { seqs_trees_templates }

      seqs_trees_templates
        .join(libs_ch, remainder: true)
        .set { seqs_trees_templates_libs }

      //seqs_trees_templates_libs.view()
      TCOFFEE_CI(seqs_trees_templates_libs, refs_ch, tc_mode)
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
