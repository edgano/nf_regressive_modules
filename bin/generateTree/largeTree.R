    #!/usr/bin/env Rscript
library('seqinr')
library('ape')

args <- commandArgs(trailingOnly = TRUE)

fastaFile <- seqinr::read.fasta(file = args[1])

    #how many fasta sequences
length(fastaFile)
    #get the seqs IDs
seqNames<-seqinr::getName(fastaFile)
    #define numbers of leave
n <- length(fastaFile)
    #produce the tree
rt<-ape::rmtree(1, n, rooted = TRUE, tip.label = seqNames, br = runif)
    #convert & save to newick
ape::write.tree(rt, file = "out.dnd", append = FALSE,digits = 2, tree.names = FALSE)
