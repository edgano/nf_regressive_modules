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
include TCOFFEE_DEFAULT         from './modules/tcoffeeModes.nf'   
include TCOFFEE_QUICKALN        from './modules/tcoffeeModes.nf'   
include TCOFFEE_MCOFFEE         from './modules/tcoffeeModes.nf'   
include TCOFFEE_ACCURATE        from './modules/tcoffeeModes.nf'   
include TCOFFEE_FMCOFFEE        from './modules/tcoffeeModes.nf'   
include TCOFFEE_PSICOFFEE       from './modules/tcoffeeModes.nf'   
include TCOFFEE_EXPRESSO        from './modules/tcoffeeModes.nf'   
include TCOFFEE_PROCOFFEE       from './modules/tcoffeeModes.nf'   
include TCOFFEE_3DCOFFEE        from './modules/tcoffeeModes.nf'   
include TCOFFEE_TRMSD           from './modules/tcoffeeModes.nf'   
include TCOFFEE_RCOFFEE         from './modules/tcoffeeModes.nf'   
include TCOFFEE_RCOFFEE_CONSAN  from './modules/tcoffeeModes.nf'   
include TCOFFEE_3DALIGN         from './modules/tcoffeeModes.nf'   
include TCOFFEE_3DMALIGN        from './modules/tcoffeeModes.nf'   

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

// tokenize params 
//pair_method = params.pairMethods.tokenize(',')

/*    main script flow    */
workflow pipeline {
    seqs_ch
      .join(templates_ch, remainder: true)
      .set { seqs_templates }

    seqs_templates
      .join(libs_ch, remainder: true)
      //.view()
      .set { seqs_templates_libs }

    TCOFFEE_DEFAULT(seqs_templates_libs,"NA")
    TCOFFEE_QUICKALN(seqs_templates_libs,"NA")
    TCOFFEE_MCOFFEE(seqs_templates_libs,"NA")
    //TCOFFEE_ACCURATE(seqs_templates_libs,"NA")
    TCOFFEE_FMCOFFEE(seqs_templates_libs,"NA")
    TCOFFEE_PSICOFFEE(seqs_templates_libs,"NA")
    TCOFFEE_EXPRESSO(seqs_templates_libs,"NA")
    TCOFFEE_PROCOFFEE(seqs_templates_libs,"NA")
    TCOFFEE_3DCOFFEE(seqs_templates_libs,"NA")
    TCOFFEE_TRMSD(seqs_templates_libs,"NA")
    //TCOFFEE_RCOFFEE(seqs_templates_libs,"NA")
    //TCOFFEE_RCOFFEE_CONSAN(seqs_templates_libs,"NA")
    TCOFFEE_3DALIGN(seqs_templates_libs,"NA")
    TCOFFEE_3DMALIGN(seqs_templates_libs,"NA")
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