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

params.trees=null
params.tree_methods="MBED"

// Input the templates
params.templates="/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/sh3.template_file"

//Input for the pdb files
params.pdb="/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/PDBs/*.pdb"

params.pairMethods = "sap_pair,tmalign_pair,mustang_pair,slow_pair"         //sap, TMalign, mustang

params.libs=""

// params for tcofee methods
params.params4tcoffee = ''

// output directory
params.outdir = "$baseDir/results_test"

          //uniref50, pdb or path
params.db = "pdb"        
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
         Tree methods                                   : ${params.tree_methods}
         Tcoffee methods                                : ${params.tc_modes}
         Template files                                 : ${params.templates}
         PDB files                                      : ${params.pdbFiles}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()

// import analysis pipelines 
include LIBRARY_GENERATION      from './modules/preprocess.nf'
include TREE_GENERATION         from './modules/treeGeneration'    
include TCOFFEE_3DALIGN         from './modules/tcoffeeModes.nf'
include TCOFFEE_3DMALIGN        from './modules/tcoffeeModes.nf'
include REG_3DALIGN             from './modules/regressiveAlignments.nf'
include REG_3DMALIGN            from './modules/regressiveAlignments.nf'
include IRMSD                   from './modules/evaluateAlignment.nf'


seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

if ( params.templates ) {
  Channel
  .fromPath(params.templates, checkIfExists: true)
  .map { item -> [ item.simpleName, item] }
  .set { templates_ch }
}else { 
  Channel
    .empty()
    .set { templates_ch }
}

if ( params.pdb ) {
  Channel
  .fromPath(params.pdb, checkIfExists: true)
  .map { item -> [ item.simpleName , item] }
  .groupTuple()
  .set { pdbFiles}
}else { 
  Channel
    .empty()
    .set { pdbFiles }
}

if ( params.pairFile ) {
  Channel
  .fromPath(params.pairFile, checkIfExists: true)
  .map { item -> [ item.simpleName , item] }
  .groupTuple()
  .set { pair_file}
}else { 
  Channel
    .empty()
    .set { pair_file }
}

if ( params.libs ) {
  Channel
  .fromPath(params.libs, checkIfExists: true)
  .map { item -> [ item.simpleName , item] }
  .groupTuple()
  .set { libs_ch }
}else { 
  Channel
    .empty()
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
bucket_list = params.buckets.toString().tokenize(',')     //int to string

/*    main script flow    */
workflow pipeline {
    seqs_ch
      .join(templates_ch, remainder: true)
      .set { seqs_templates }

    LIBRARY_GENERATION(seqs_templates)

    seqs_templates
      .join(LIBRARY_GENERATION.out.sap_lib, remainder: true)
      //.view()
      .set { seqs_templates_libs }
        //TODO -> is needed to merge all the LIBS   ???

    TCOFFEE_3DALIGN(seqs_templates_libs,"NA") //COMPUTE_3D_ALING()

    TCOFFEE_3DMALIGN(seqs_templates_libs,"NA")//COMPUTE_3DM_ALING()

    if (!params.trees){
      TREE_GENERATION (seqs_ch, tree_method) 
      seqs_ch
        .cross(TREE_GENERATION.out)
        .map { it -> [ it[1][0], it[1][1], it[1][2] ] }
        .set { seqs_and_trees }
    }else{
      seqs_ch
        .cross(trees)
        .map { it -> [ it[1][0], it[1][1], it[1][2] ] }
        .set { seqs_and_trees }
    }
    seqs_and_trees
      .join(seqs_templates_libs, remainder: true)
      .view()
      .set { seqs_trees_templates_libs }

    REG_3DALIGN(seqs_trees_templates_libs,"1000")                                  

    //REG_3DMALIGN(seqs_trees_templates_libs,"1000") 

    REG_3DALIGN.out.alignmentFile
      .join(templates_ch, remainder: true)
      .set { aln_templates }

    IRMSD(aln_templates,"3dAlign")
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