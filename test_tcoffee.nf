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
//params.seqs = "/users/cn/egarriga/datasets/homfam/combinedSeqs/{${top20fam}}.fa"
params.seqs ="${baseDir}/test/three.fa"

//params.refs = "/users/cn/egarriga/datasets/homfam/refs/{${top20fam}}.ref"
params.refs ="${baseDir}/test/three.ref"

//params.trees ="/Users/edgargarriga/CBCRG/nf_regressive_modules/results/trees/*.dnd"
params.tree="${baseDir}/test/three.MBED.dnd"

params.templates="${baseDir}/test/*_ref.template_list"

params.pdbFiles="${baseDir}/test/*.pdb"

//default,quickaln,mcoffee,fmcoffee,accurate,psicoffee,expresso,procoffee,3dcoffee,trmsd,rcoffee
    //accurate,expresso         Impossible to find EXPRESSO Templates Check that your blast server is properly installed [See documentation][FATAL:T-COFFEE] 
    //mcoffee -> poa is needed
    //trmsd --> templates
    //rcoffee --> RNAplfold
params.tc_modes = "default,quickaln,fmcoffee,psicoffee,procoffee,3dcoffee"
//default,quickaln,fmcoffee,psicoffee,procoffee,3dcoffee"

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
include TREE_GENERATION from './modules/treeGeneration'        params(params)
include TCOFFEE_ANALYSIS from './modules/prog_analysis'    params(params)

// Channels containing sequences
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

if ( params.templates ) {
  Channel
  .fromPath(params.templates)
//  .map { item -> [ item.baseName.replace('_ref','') , item] }
  .set { templates }
}else { 
  Channel.empty().set { templates }
}

if ( params.pdbFiles ) {
  Channel
  .fromPath(params.pdbFiles)
//  .map { item -> [ item.simpleName , item] }
  .collect()
  .set { pdbFiles }
}else { 
  Channel.empty().set { pdbFiles }
}

// tokenize params 
tcoffee_mode = params.tc_modes.tokenize(',')

/*    main script flow    */
workflow pipeline {
      TCOFFEE_ANALYSIS(seqs_ch, tcoffee_mode)
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
