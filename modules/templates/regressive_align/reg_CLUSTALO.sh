valgrind t_coffee -reg -reg_method clustalo_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy -n_core=1 \
         -outfile ${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln
