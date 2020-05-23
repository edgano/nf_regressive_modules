export NO_MAFFT_BINARIES=1

t_coffee -other_pg seq_reformat -in ${guide_tree} -input newick -in2 ${seqs} -input2 fasta_seq -action +newick2mafftnewick >> ${id}.mafftnewick

newick2mafft.rb 1.0 ${id}.mafftnewick > ${id}.mafftbinary

mafft --retree 1 --anysymbol --treein ${id}.mafftbinary ${seqs} > ${id}.prog.${align_method}.with.${tree_method}.tree.aln
