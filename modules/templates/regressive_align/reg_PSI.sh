export blast_server_4_CLTCOFFEE=LOCAL
export protein_db_4_CLTCOFFEE=${params.database_path}
export VERBOSE_4_DYNAMIC=1

t_coffee -reg -reg_method psicoffee_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -cache ${params.blastOutdir} \
         -outfile ${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln