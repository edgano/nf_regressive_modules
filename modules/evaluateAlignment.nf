#!/bin/bash nextflow
params.outdir = 'results'

process EVAL_ALIGNMENT {
    container 'edgano/tcoffee:slave'
    tag "EVAL_ALIGNMENT on $id"
    publishDir "${params.outdir}/score"

    input:
    val align_type
    val id
    tuple file (test_alignment), val (id), file (ref_alignment)
    val align_method
    val tree_method
    val bucket_size

    output:
    file("*.sp")
    path "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc", emit: tcScore
    file("*.col") 

    script:
    """
    ## Sum-of-Pairs Score ##
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.sp"

    ## Total Column Score ##	
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode tc \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc"

    ## Column Score ##
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode column \
            | grep -v "seq1" | grep -v '*' | \
              awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.col"       
    """
}

process EASEL_INFO {
    container 'edgano/hmmer:latest'
    tag "EASEL_INFO on $id"
    publishDir "${params.outdir}/easel"

    input:
    val align_type
    val id
    tuple file (test_alignment), val (id), file (ref_alignment)
    val align_method
    val tree_method
    val bucket_size

    output:
    file("*.easel_INFO")
    file("*.avgLen")
    file("*.avgId")
      
     shell:
     '''
     esl-alistat !{test_alignment} > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_INFO
     awk -F : '{ if (\$1=="Average length") print \$2}' !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_INFO | sed 's/ //g' > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.avgLen 
     awk -F : '{ if (\$1=="Average identity") print substr(\$2, 1, length(\$2)-1)}' !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_INFO | sed 's/ //g' > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.avgId 

     ## awk 'NR > 8 && $1 !~/\\// { sum+= $3 } END {print "SUM: "sum"\\nAVG: "sum/(NR-9)}' !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_INFO > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_AVG
     ## the first && is to skip first lines and the last one. The AVG is done -8 all the time execpt for the END print to "erase" the last "//" too.
     '''
}

process HOMOPLASY {
    container 'edgano/base:latest'
    tag "HOMOPLASY on $id"
    publishDir "${params.outdir}/homoplasy"

    input:
    val align_type
    val id
    tuple file (test_alignment), val (id), file (ref_alignment)
    val align_method
    val tree_method
    val bucket_size
    file (homoplasy)        

    output:
    file("*.homo")
    file("*.w_homo")
    file("*.w_homo2")
    file("*.len")
    file("*.ngap")
    file("*.ngap2")

    script:
  """    
    ## homo
    awk -F : '{ if (\$1=="HOMOPLASY") print \$2}' ${homoplasy} > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.homo
    ## w_homo
    awk -F : '{ if (\$1=="WEIGHTED_HOMOPLASY") print \$2}' ${id}.homoplasy > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo
    ## w_homo2
    awk -F : '{ if (\$1=="WEIGHTED_HOMOPLASY2") print \$2}' ${id}.homoplasy > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo2
    ## len
    awk -F : '{ if (\$1=="LEN") print \$2}' ${id}.homoplasy > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.len
    ## ngap
    awk -F : '{ if (\$1=="NGAP") print \$2}' ${id}.homoplasy > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.ngap
    ## ngap2
    awk -F : '{ if (\$1=="NGAP2") print \$2}' ${id}.homoplasy > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.ngap2 
    """
}

process METRICS {
    container 'edgano/base:latest'
    tag "METRICS on $id"
    publishDir "${params.outdir}/metrics"

    input:
        val align_type
        val id
        tuple file (test_alignment), val (id), file (ref_alignment)
        val align_method
        val tree_method
        val bucket_size
        file (metricsFile)

    output:
      file("*.metrics")
      file("*.realtime")
      file("*.rss")
      file("*.peakRss")
      file("*.vmem")
      file("*.peakVmem")
      file("*.metrics") 

    script:
    """    
        ## realtime > Task execution time i.e. delta between completion and start timestamp i.e. compute wall-time
    awk -F = '{ if (\$1=="realtime") print \$2}' ${metricsFile} > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.realtime
        ## rss > Real memory (resident set) size of the process
    awk -F = '{ if (\$1=="rss") print \$2}' ${metricsFile}> ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.rss
        ## peakRss > Peak of real memory
    awk -F = '{ if (\$1=="peak_rss") print \$2}' ${metricsFile} > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.peakRss
        ## vmem > Virtual memory size of the process
    awk -F = '{ if (\$1=="vmem") print \$2}' ${metricsFile} > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.vmem
        ## peakVmem > Peak of virtual memory
    awk -F = '{ if (\$1=="peak_vmem") print \$2}' ${metricsFile} > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.peakVmem
    
    mv ${metricsFile} ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.metrics
    """
}

process GAPS_PROGRESSIVE {
    container 'edgano/base:latest'
    tag "GAPS_PROG on $id"
    publishDir "${params.outdir}/metrics"

    input:
    val align_type
    val id
    file test_alignment
    val align_method
    val tree_method
    val bucket_size

    output:
        file("*.totGap")
        file("*.numSeq")
        file("*.alnLen")

    script:
    """    
#!/usr/bin/env python
from Bio import SeqIO
from decimal import *
import os
gap = '-'
globalGap = 0
avgGap = 0
auxGap = 0   
totGapName= "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.totGap"
numbSeqName= "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.numSeq"
alnLenName= "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.alnLen"
totGapFile= open(totGapName,"w+")
numSeqFile= open(numbSeqName,"w+")
alnLenFile= open(alnLenName,"w+")
record = list(SeqIO.parse("${test_alignment}", "fasta"))
for sequence in record:
    ## print(sequence.seq)
    auxGap = sequence.seq.count(gap)
    globalGap += auxGap
avgGap = Decimal(globalGap) / Decimal(len(record))
print "NumSeq: ",len(record)," GlobalGap: ",globalGap," AVG_Gap:",round(avgGap,3)
totGapFile.write(str(globalGap))
alnLenFile.write(str(len(record[0])))
numSeqFile.write(str(len(record)))
totGapFile.close()
alnLenFile.close()
numSeqFile.close()
    """
}
