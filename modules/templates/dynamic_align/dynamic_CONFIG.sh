t_coffee -reg -reg_method dynamic_msa \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -dynamic ${dynamic_size} \
         -reg_homoplasy \
         -dynamic_config $PWD/config.txt \
         -blast_server LOCAL \
         -protein_db=${params.db} \
         -outfile ${id}.dynamic.${bucket_size}.dynamicSize.${dynamic_size}.${align_method}.with.${tree_method}.tree.aln