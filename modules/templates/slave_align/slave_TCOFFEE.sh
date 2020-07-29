t_coffee -reg -reg_method t_coffee_msa \
         -reg_tree ${guide_tree} \
         -child_tree ${slave_method} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.slave_${bucket_size}.${align_method}.with.${tree_method}_${slave_method}.tree.aln
