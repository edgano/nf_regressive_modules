export NO_MAFFT_BINARIES=1
export MAFFT_BINARIES=''

t_coffee -reg -reg_method mafft_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln
