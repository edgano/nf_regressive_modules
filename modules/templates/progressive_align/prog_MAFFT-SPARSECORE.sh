replace_U.pl ${seqs}

mafft-sparsecore.rb -i ${seqs} > ${id}.prog.${align_method}.with.NO_TREE.tree.aln