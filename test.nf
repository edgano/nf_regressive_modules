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
params.seqs = "/users/cn/egarriga/datasets/homfam/combinedSeqs/*.fa"

//params.trees ="/Users/edgargarriga/CBCRG/nf_regressive_modules/results/trees/*.dnd"
params.trees = false
                      //MAFFT-DPPARTTREE0,FAMSA-SLINK,MBED,MAFFT-PARTTREE0
params.tree_methods = "FAMSA-SLINK"

// output directory
params.outdir = "$baseDir/results"


log.info """\
         PIPELINE  ~  version 0.1"
         ======================================="
         Input trees (NEWICK)                           : ${params.trees}
         Tree methods                                   : ${params.tree_methods}
         --##--
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()

// import analysis pipelines
include TREE_GENERATION from './modules/treeGeneration'   params(params)
include SACKIN_INDEX from './modules/evaluateAlignment'            params(params)

seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }
// Channels for user provided trees or empty channel if trees are to be generated [OPTIONAL]
if ( params.trees ) {
  trees = Channel.fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
}else { 
  Channel.empty().set { trees }
}

// tokenize params 
tree_method = params.tree_methods.tokenize(',')

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

    
      SACKIN_INDEX(seqs_and_trees)
    
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
