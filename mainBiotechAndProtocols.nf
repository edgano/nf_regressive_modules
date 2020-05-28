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
//params.seqs ="/users/cn/egarriga/datasets/homfam/combinedSeqs/{${seq2improve}}.fa"
top20fam="gluts,myb_DNA-binding,tRNA-synt_2b,biotin_lipoyl,hom,ghf13,aldosered,hla,Rhodanese,PDZ,blmb,rhv,p450,adh,aat,rrm,Acetyltransf,sdr,zf-CCHH,rvp"

// input sequences to align in fasta format
params.seqs = "/users/cn/egarriga/datasets/homfam/combinedSeqs/{${top20fam}}.fa"

params.refs = "/users/cn/egarriga/datasets/homfam/refs/*.ref"

//params.trees ="/Users/edgargarriga/CBCRG/nf_regressive_modules/results/trees/*.dnd"
params.trees = false
                      //CLUSTALO,FAMSA,MAFFT-FFTNS1
params.align_methods = "CLUSTALO,FAMSA,MAFFT-FFTNS1" 
                      //MAFFT-DPPARTTREE0,FAMSA-SLINK,MBED,MAFFT-PARTTREE
params.tree_methods = "MAFFT-PARTTREE" 

params.buckets = "1000"


params.biotech = false
params.reg_protocols = false
params.tcoffee_protocols = false

params.evaluate=true
params.homoplasy=false
params.gapCount=false
params.metrics=false
params.easel=false

// output directory
params.outdir = "$baseDir/resultsBiotechAndProtocols"

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
         Generate Nature Biotech analysis               : ${params.biotech}
         Generate Regressive Protocols analysis         : ${params.reg_protocols}
         Generate TCoffee Protocols analysis            : ${params.tcoffee_protocols}
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
include TREE_GENERATION             from './modules/treeGeneration'         params(params)
include NAT_BIOTECH_ANALYSIS        from './modules/papers_analysis'        params(params)
include REG_PROTOCOLS_ANALYSIS      from './modules/papers_analysis'        params(params)
include TCOFFEE_PROTOCOLS_ANALYSIS  from './modules/papers_analysis'        params(params)

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

// tokenize params 
tree_method = params.tree_methods.tokenize(',')
align_method = params.align_methods.tokenize(',')
bucket_list = params.buckets.toString().tokenize(',')       //int to string

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

    if (params.biotech){
      //TODO define aligners & tree
        // CO,FFTNS1,UPP,SPARSECORE,GINSI
        // PARTREE0,MBED
        // 1000
        // Homfam dataset
      NAT_BIOTECH_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)
    }
    if (params.reg_protocols){
      //TODO define aligners & tree
        // CO,FFTNS1,GINSI
        // PARTREE0,MBED
        // 1000
        // gluts.fasta
        
        //t_coffee -reg -seq gluts.fasta -reg_nseq 1000 -reg_tree mbed -reg_method clustalo_msa -outfile gluts.aln -outtree gluts.mbed
        //t_coffee -reg -seq gluts.fasta -reg_nseq 1000 -reg_tree mbed -reg_method mafftginsi_msa -outfile gluts.aln -outtree gluts.mbed
        //t_coffee -reg -seq gluts.fasta -reg_nseq 1000 -reg_tree parttree -reg_method mafftfftnsi_msa -outfile gluts.aln -outtree gluts.parttree
      REG_PROTOCOLS_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)
    }
    if (params.tcoffee_protocols){
      //TODO define parameters

      TCOFFEE_PROTOCOLS_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)
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
