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
params.seqs = ""  //  "/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/sh3.fasta"

params.refs = ""  //  "/users/cn/egarriga/datasets/homfam/refs/*.ref"

//params.trees ="/users/cn/egarriga/datasets/homfam/trees/*.{FAMSA,CLUSTALO,MAFFT_PARTTREE}.dnd"
params.trees = false

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
params.tree_methods = "FAMSA-SLINK"      

//blast call cached
params.generateBlast=false
params.blastOutdir="$baseDir/blast"

// ## TCOFFEE
                  //3DALIGN,3DCOFFEE,3DMALIGN,ACCURATE,DEFAULT,EXPRESSO,FMCOFFEE,MCOFFEE,PROCOFFEE,PSICOFFEE,QUICKALN,RCOFFEE_CONSAN,RCOFFEE,TRMSD"
  //need template -> 3DMALIGN
params.tc_modes = "DEFAULT"
params.templates = ''//   '/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/sh3.template_file'
params.pdb = ''     //   '/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/PDBs/*.pdb'
params.libs = ''
params.pairFile = ''
params.params4tcoffee = ''   
params.cache_path = "${params.blastOutdir}"

          //uniref50, pdb or path
params.db = "pdb"        
          // define default database path
uniref_path = "/users/cn/egarriga/datasets/db/uniref50.fasta"   // cluster path
pdb_path = "/database/pdb/pdb_seqres.txt"                       // docker path

  
params.tcoffee_align = true            

// output directory
params.outdir = "$baseDir/results"

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
         --##--
         Generate TCoffee alignments                    : ${params.tcoffee_align}
                  TCoffee modes                         : ${params.tc_modes}
         --##--
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()



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
<<<<<<< HEAD
//tcoffee_mode = params.tc_modes.tokenize(',')
//pair_method = params.pairMethods.tokenize(',')
tree_method = params.tree_methods.tokenize(',' as char)
bucket_list = params.buckets.toString().tokenize(',')     //int to string
=======
tree_method = params.tree_methods.tokenize(',')
tc_mode = params.tc_modes.tokenize(',')
>>>>>>> bebf4fa34f11a60205ccf1ebaf050708dfebd9de


 // import analysis pipelines
include TREE_GENERATION from './modules/treeGeneration'        params(params)
include TCOFFEE_CI from './modules/prog_analysis'      params(params)

/*   main script flow     */
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

    if (params.tcoffee_align){

      seqs_and_trees
        .join(templates_ch, remainder: true)
        .set { seqs_trees_templates }

      seqs_trees_templates
        .join(libs_ch, remainder: true)
        .set { seqs_trees_templates_libs }

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
