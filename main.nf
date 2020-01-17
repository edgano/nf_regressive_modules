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

 */

/*
 * defaults parameter definitions
 */

// input sequences to align in fasta format
params.seqs = "$baseDir/data/"

// output directory
params.output = "$baseDir/results"


log.info """\
         PIPELINE  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Output directory (DIRECTORY)                   : ${params.output}
         """
         .stripIndent()


// Channels containing sequences
if ( params.seqs ) {
  Channel
  .fromPath(params.seqs)
  .map { item -> [ item.baseName, item] }
  .into { seqs }
}

process getAlnGaps {
    tag "${id}"
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
      set val(id), file(aln) from seqs

    output:
     set val(id), file("*") into gapsOut

    script:
    """
    """
}