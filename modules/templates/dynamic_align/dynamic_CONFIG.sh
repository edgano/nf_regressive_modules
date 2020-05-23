export blast_server_4_CLTCOFFEE=LOCAL
export protein_db_4_CLTCOFFEE=${params.database_path}
export NO_MAFFT_BINARIES=1

t_coffee -reg -reg_method dynamic_msa \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -dynamic ${dynamicX} \
         -reg_homoplasy \
         -dynamic_config $PWD/config.txt \
         -outfile ${id}.dynamic_${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.with.${tree_method}.tree.aln
