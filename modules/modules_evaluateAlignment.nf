#!/bin/bash nextflow
params.outdir = 'results'

process EVAL_ALIGNMENT {
    container 'edgano/tcoffee:pdb'
    tag "EVAL_ALIGNMENT on $id"
    publishDir "${params.outdir}/score/tc", mode: 'copy', overwrite: true, pattern: '*.tc'
    publishDir "${params.outdir}/score/sp", mode: 'copy', overwrite: true, pattern: '*.sp'
    publishDir "${params.outdir}/score/col", mode: 'copy', overwrite: true, pattern: '*.col'

    input:
    val align_type
    tuple  val (id), file (test_alignment), file (ref_alignment)
    val align_method
    val tree_method
    val bucket_size

    output:
    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.sp"), emit: spScore

    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc"), emit: tcScore

    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.col"), emit: colScore
   
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
    publishDir "${params.outdir}/easel", mode: 'copy', overwrite: true

    input:
    val align_type
    tuple  val (id), file (test_alignment)
    val align_method
    val tree_method
    val bucket_size

    output:
    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path("*.easel_INFO"), \
    path("*.avgLen"), \
    path("*.avgId"), emit: easelFiles
      
     shell:
     '''
     esl-alistat !{test_alignment} > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_INFO
     awk -F : '{ if (\$1=="Average length") printf "%s", \$2}' !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_INFO | sed 's/ //g' > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.avgLen 
     awk -F : '{ if (\$1=="Average identity") printf "%s", substr(\$2, 1, length(\$2)-1)}' !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_INFO | sed 's/ //g' > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.avgId 

     ## awk 'NR > 8 && $1 !~/\\// { sum+= $3 } END {print "SUM: "sum"\\nAVG: "sum/(NR-9)}' !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_INFO > !{id}.!{align_type}.!{bucket_size}.!{align_method}.with.!{tree_method}.tree.easel_AVG
     ## the first && is to skip first lines and the last one. The AVG is done -8 all the time execpt for the END print to "erase" the last "//" too.
     '''
}

process HOMOPLASY {
    container 'edgano/base:latest'
    tag "HOMOPLASY on $id"
    publishDir "${params.outdir}/homoplasy", mode: 'copy', overwrite: true

    input:
    val align_type
    tuple  val (id), file (test_alignment)
    val align_method
    val tree_method
    val bucket_size
    file homoplasy     

    output:
    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("*.homo"), \
    path("*.w_homo"), \
    path("*.w_homo2"), \
    path("*.len"), \
    path("*.ngap"), \
    path("*.ngap2"), emit: homoFiles

    script:
  """  
        ## remove whitespace  
    cat ${homoplasy} | tr -d ' ' > aux.txt 
        ## homo     
    awk -F : '{ if (\$1=="HOMOPLASY") printf "%s", \$2}' aux.txt  > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.homo 
        ## w_homo
    awk -F : '{ if (\$1=="WEIGHTED_HOMOPLASY") printf "%s", \$2}' aux.txt  > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo
        ## w_homo2
    awk -F : '{ if (\$1=="WEIGHTED_HOMOPLASY2") printf "%s", \$2}' aux.txt > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo2
        ## len
    awk -F : '{ if (\$1=="LEN") printf "%s", \$2}' aux.txt > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.len
        ## ngap
    awk -F : '{ if (\$1=="NGAP") printf "%s", \$2}' aux.txt > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.ngap
        ## ngap2
    awk -F : '{ if (\$1=="NGAP2") printf "%s", \$2}' aux.txt > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.ngap2 
    """
}

process METRICS {
    container 'edgano/base:latest'
    tag "METRICS on $id"
    publishDir "${params.outdir}/metrics", mode: 'copy', overwrite: true

    input:
    val align_type
    tuple  val (id), file (test_alignment)
    val align_method
    val tree_method
    val bucket_size
    file metricsFile

    output:
    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path("*.realtime"), \
    path("*.rss"), \
    path("*.peakRss"), \
    path("*.vmem"), \
    path("*.peakVmem"), emit: metricFiles

    script:
    """    
        ## remove whitespace  
    cat ${metricsFile} | tr -d ' ' > aux.txt 
        ## realtime > Task execution time i.e. delta between completion and start timestamp i.e. compute wall-time
    awk -F = '{ if (\$1=="realtime") printf "%s", \$2}' aux.txt  > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.realtime
        ## rss > Real memory (resident set) size of the process
    awk -F = '{ if (\$1=="rss") printf "%s", \$2}' aux.txt > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.rss
        ## peakRss > Peak of real memory
    awk -F = '{ if (\$1=="peak_rss") printf "%s", \$2}' aux.txt  > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.peakRss
        ## vmem > Virtual memory size of the process
    awk -F = '{ if (\$1=="vmem") printf "%s", \$2}' aux.txt  > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.vmem
        ## peakVmem > Peak of virtual memory
    awk -F = '{ if (\$1=="peak_vmem") printf "%s", \$2}' aux.txt  > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.peakVmem
    
    ## mv ${metricsFile} ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.metrics
    """
}

process GAPS_PROGRESSIVE {
    container 'edgano/base:latest'
    tag "GAPS_PROG on $id"
    publishDir "${params.outdir}/gaps", mode: 'copy', overwrite: true

    input:
    val align_type
    tuple  val (id), file (test_alignment)
    val align_method
    val tree_method
    val bucket_size

    output:
    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path("*.totGap"), \
    path("*.numSeq"), \
    path("*.alnLen"), emit: gapFiles

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

process IRMSD{
    container 'edgano/tcoffee:protocols'
    tag "$align_type on $id "
    publishDir "${params.outdir}/iRMSD/${id}"
    
    input:
        val align_type
        set val(id), path(alignment), path(template)
        val align_method
        val tree_method
        val bucket_size
    
    output:
        path "*.irmsd", emit: irmsd
    
    script:
    """
    t_coffee -other_pg irmsd ${alignment} -template_file ${template} -io_format s > ${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.irmsd
    """
}

process SACKIN_INDEX {
    container 'edgano/r_base:latest'
    tag " on $id - $tree_method"
    publishDir "${params.outdir}/sackin"
    
    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    
    output:
    path("*_sackin")

    script:
    """
    #!/usr/bin/env Rscript

    library(apTreeshape)

    tree <- as.treeshape(read.tree("${guide_tree}"))

            ## Index of Sackin of a PDA tree :    info https://rdrr.io/cran/apTreeshape/man/shape.statistic.html
    sackin_pda <- sackin(tree,norm="pda")
    sackin_yule <- sackin(tree,norm="yule")

    sackin.test(tree,model='pda')
    sackin.test(tree,model='yule')

    write(sackin_pda, file = "${id}_${tree_method}.pda_sackin", append = FALSE, sep = " ")
    """
}
