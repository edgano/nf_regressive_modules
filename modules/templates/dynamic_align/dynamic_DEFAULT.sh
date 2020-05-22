export blast_server_4_CLTCOFFEE=LOCAL
export protein_db_4_CLTCOFFEE=${params.database_path}

t_coffee -reg -reg_method dynamic_msa \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -dynamic ${dynamicX} \
         -reg_homoplasy \
         -outfile ${id}.dynamic_${bucket_size}.dynamicX.${dynamicX}.DEFAULT.with.${tree_method}.tree.aln