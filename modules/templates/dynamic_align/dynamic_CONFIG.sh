  export blast_server_4_CLTCOFFEE=LOCAL
  export protein_db_4_CLTCOFFEE=${params.db}

t_coffee -reg -reg_method dynamic_msa \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -dynamic ${dynamic_size} \
         -reg_homoplasy \
         -dynamic_config $PWD/config.txt \
         -outfile ${id}.dynamic.${bucket_size}.dynamicSize.${dynamic_size}.${align_method}.with.${tree_method}.tree.aln