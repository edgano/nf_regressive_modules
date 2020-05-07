t_coffee -reg -reg_method clustalo_msa \
     -seq ${seqs} \
     -reg_tree ${guide_tree} \
     -child_tree ${slave_method} \
     -reg_nseq ${bucket_size} \
     -reg_homoplasy \
     -outfile ${id}.slave.${bucket_size}.${align_method}.with.${tree_method}.tree.slave.${slave_method}.aln
              