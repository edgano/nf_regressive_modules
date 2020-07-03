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

// input sequences to align in fasta format

// input the alignments file
params.seqs="/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/sh3.fasta"

// Input the templates
params.templates="/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/sh3.template_file"

//Input for the pdb files
params.pdb="/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/PDBs/*.pdb"

//##--##--#   INTRAMOL_MATRIX_GENERATION    #--##--##
  params.threedTreeMode       = "4"     // 1-6
  params.threedTreeNoWeights  = "2"     // (2=square; 3=cubic)
  params.threedTreeModeExp    = "1"

  params.replicatesNum        = "1"
  params.evaluate3DVal        = "contacts"    //"contacts,distances,strike"
//##--##--##--##--####--##--##--##--####--##--##--##--##

params.pairMethods = "sap_pair,tmalign_pair,mustang_pair,slow_pair"         //sap, TMalign, mustang

params.pairFile=""

params.libs=""

//default,quickaln,mcoffee,fmcoffee,psicoffee,expresso,procoffee,3dcoffee,trmsd,rcoffee,rcoffee_consan"
    //accurate    --> TODO   
    //3d_align, 3dM_align    
params.tc_modes = "default"//,quickaln,mcoffee,fmcoffee,psicoffee,expresso,procoffee"//3dcoffee,trmsd"

// params for tcofee methods
params.params4tcoffee = ''

// output directory
params.outdir = "$baseDir/results_test"

          //uniref50, pdb or path
params.db = "pdb"        
          // define database path
uniref_path = "/users/cn/egarriga/datasets/db/uniref50.fasta"   // cluster path
pdb_path = "/database/pdb/pdb_seqres.txt"                       // docker path

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
         Tcoffee methods                                : ${params.tc_modes}
         Template files                                 : ${params.templates}
         PDB files                                      : ${params.pdbFiles}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()

// import analysis pipelines
include CLANS_ANALYSIS from './modules/leila_analysis'       params(params) 
include TCOFFEE_TEST from './modules/prog_analysis'         params(params)
include TCOFFEE_ANALYSIS from './modules/prog_analysis'     params(params) 

seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

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
  .set { libs }
}else { 
  Channel.empty()
    .set { libs }
}

if ( params.trees ) {
  trees_ch = Channel.fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
}else { 
  Channel.empty().set { trees_ch }
}
// tokenize params 
//tcoffee_mode = params.tc_modes.tokenize(',')
//pair_method = params.pairMethods.tokenize(',')
tree_method = params.tree_methods.tokenize(',' as char)
bucket_list = params.buckets.toString().tokenize(',')     //int to string

/*    main script flow    */
workflow pipeline {

      //TCOFFEE_ANALYSIS(seqs_ch, templates, pdbFiles, pair_file, libs, tcoffee_mode, pair_method)
      //TCOFFEE_TEST(seqs_ch, templates, pdbFiles, pair_file, libs, tcoffee_mode, pair_method)
      //CLANS_ANALYSIS(seqs_ch, trees_ch, tree_method, bucket_list,templates_ch)

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