replace_U.pl ${seqs}

run_upp.py -s ${seqs} \
           -m amino \
           -x 1 \
           -o ${id}.prog.${align_method}

mv ${id}.prog.${align_method}_alignment.fasta ${id}.prog.${align_method}.with.NO_TREE.tree.aln