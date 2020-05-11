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

params.align_method = "CLUSTALO"

params.tree_method = "CLUSTALO"

params.buckets = '1000'

params.trees = true
params.regressive_align = true

// output directory
params.outdir = "$baseDir/results"

log.info """\
         PIPELINE  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (FASTA)                       : ${params.refs}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()

// import modules
include REG_ANALYSIS from './modules/reg_analysis'    params(params)

// Channels containing sequences
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

if ( params.refs ) {
  refs_ch = Channel.fromPath( params.refs ).map { item -> [ item.baseName, item] }
}

tree_methods = params.tree_method
align_methods = params.align_method
bucket_list = params.buckets

/* 
 * main script flow
 */
workflow pipeline {
    REG_ANALYSIS(seqs_ch, refs_ch, align_methods, tree_methods, bucket_list)
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
