export NO_MAFFT_BINARIES=1
export MAFFT_BINARIES=''

replace_U.pl ${seqs} 
t_coffee -reg -reg_method mafftsparsecore_msa \
        -reg_tree ${guide_tree} \
        -seq ${seqs} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        -outfile ${id}.pool_${bucket_size}.${align_method}.with.${tree_method}.tree.aln