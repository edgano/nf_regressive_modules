export NO_MAFFT_BINARIES=1
export MAFFT_BINARIES=''

t_coffee -reg -reg_method clustalo_msa \
     -seq ${seqs} \
     -reg_tree ${guide_tree} \
     -child_tree ${slave_method} \
     -reg_nseq ${bucket_size} \
     -reg_homoplasy \
     -outfile ${id}.slave_${bucket_size}.${align_method}.with.${tree_method}_${slave_method}.tree.aln
              