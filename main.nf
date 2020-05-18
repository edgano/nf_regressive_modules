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
//params.seqs = "$baseDir/data/*.fa"
params.seqs = 'https://raw.githubusercontent.com/edgano/datasets-test/homfam/seatoxin.fa' //#TODO

//params.refs = "$baseDir/data/*.ref"
params.refs = 'https://raw.githubusercontent.com/edgano/datasets-test/homfam/seatoxin.ref' //#TODO

params.trees = false
//params.trees ="/Users/edgargarriga/CBCRG/nf_regressive_modules/results/trees/seatoxin.MBED.dnd"

params.align_methods = "CLUSTALO,MAFFT-FFTNS1"

params.tree_methods = "MBED,FAMSA-SLINK"

params.buckets = "1000,2000"

params.progressive_align = false
params.regressive_align = true
params.slave_align=false
params.slave_tree_methods="mbed" //need to be lowercase -> direct to tcoffee CommandLine
params.dynamic_align=false
params.pool_align=false

params.evaluate=true
params.homoplasy=true
params.gapCount=true
params.metrics=true

// output directory
params.outdir = "$baseDir/results"

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
                   Slave Tree methods                   : ${params.slave_tree_methods}
         Generate Dynamic alignments                    : ${params.dynamic_align}
         Generate Pool alignments                       : ${params.pool_align}
         --##--
         Perform evaluation? Requires reference         : ${params.evaluate}
         Check homoplasy? Only for regressive           : ${params.homoplasy}
         Check gapCount? For progressive                : ${params.gapCount}
         Check metrics?                                 : ${params.metrics}
         --##--
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()

// import analysis pipelines
include REG_ANALYSIS from './modules/reg_analysis'    params(params)
include PROG_ANALYSIS from './modules/prog_analysis'    params(params)
include POOL_ANALYSIS from './modules/reg_analysis'    params(params)
include SLAVE_ANALYSIS from './modules/reg_analysis'    params(params)

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

tree_method = params.tree_methods.tokenize(',')
align_method = params.align_methods.tokenize(',')
bucket_list = params.buckets.tokenize(',')
slave_method = params.slave_tree_methods.tokenize(',')

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
    if (params.pool_align){
      POOL_ANALYSIS(seqs_ch, refs_ch, align_method, tree_method, bucket_list, trees)
    }
}

workflow {
  pipeline()
  
  //publish:
    //REG_ANALYSIS.out.tcScore_csv to: "_tc.csv"
}

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
  //TODO script to generate CSV from individual files
}
