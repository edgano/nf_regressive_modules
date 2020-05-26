export NO_MAFFT_BINARIES=1

t_coffee -reg -reg_method msaprobs_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln