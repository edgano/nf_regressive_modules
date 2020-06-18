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
//params.seqs ="${baseDir}/test/RNA_set.fa"   //-->for rcoffee, rcoffee_consan
params.seqs ="${baseDir}/test/three.fa"

// input the alignments file
params.msafile="/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/sh3.fasta"

// Input the templates
params.templates="/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/sh3.template_file"

//Input for the pdb files
params.pdb="/Users/edgargarriga/CBCRG/NatureProtocolDataset/Proteins/PDBs/*.pdb"

params.pairs ="/Users/edgargarriga/CBCRG/nf_regressive_modules/test/*.pair"

//##--##--#   INTRAMOL_MATRIX_GENERATION    #--##--##
  params.threedTreeMode       = "4"     // 1-6
  params.threedTreeNoWeights  = "2"     // (2=square; 3=cubic)
  params.threedTreeModeExp    = "1"

  params.replicatesNum        = "1"
  params.evaluate3DVal        = "contacts"    //"contacts,distances,strike"
//##--##--##--##--####--##--##--##--####--##--##--##--##

params.pairMethods = "sap_pair,tmalign_pair,mustang_pair,slow_pair"

//sap, TMalign, mustang

//default,quickaln,mcoffee,fmcoffee,psicoffee,expresso,procoffee,3dcoffee,trmsd,rcoffee,rcoffee_consan"
    //accurate    --> TODO       
params.tc_modes = "accurate"//"default,quickaln,mcoffee,fmcoffee,psicoffee,expresso,procoffee,3dcoffee,trmsd"


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
include CLANS_ANALYSIS from './modules/prog_analysis'    params(params) 
include TCOFFEE_TEST from './modules/prog_analysis'    params(params)
include TCOFFEE_ANALYSIS from './modules/prog_analysis'    params(params)

seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

if ( params.msafile ) {
  Channel
  .fromPath(params.msafile)
  .map { item -> [ item.simpleName , item] }
  .set { fastaAln }
}

if ( params.templates ) {
  Channel
  .fromPath(params.templates)
  .map { item -> [ item.simpleName, item] }
  .set { templates}
}

if ( params.pdb ) {
  Channel
  .fromPath(params.pdb)
  .map { item -> [ item.simpleName , item] }
  .groupTuple()
  .set { pdbFiles}
}

if ( params.pairs ) {
  Channel
  .fromPath(params.pairs)
  .map { item -> [ item.simpleName , item] }
  .view()
  .set { pairFile}
}

fastaAln
  .combine(templates, by: 0)
  .set{ aln_templates}

/*fastaAln
  .combine(templates, by: 0)
  .combine(pdbFiles, by: 0)
  .set{ aln_templates_pdb}
*/

// tokenize params 
tcoffee_mode = params.tc_modes.tokenize(',')
pair_method = params.pairMethods.tokenize(',')

/*    main script flow    */
workflow pipeline {
      TCOFFEE_ANALYSIS(seqs_ch, tcoffee_mode, aln_templates, pair_method)
      //TCOFFEE_TEST(seqs_ch, tcoffee_mode, aln_templates, pair_method, pairFile)
      //CLANS_ANALYSIS(seqs_ch, tcoffee_mode, aln_templates, pair_method)
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
