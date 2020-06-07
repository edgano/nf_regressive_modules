export blast_server_4_CLTCOFFEE=LOCAL
export protein_db_4_CLTCOFFEE=${params.database_path}
export NO_MAFFT_BINARIES=1
export VERBOSE_4_DYNAMIC=1

t_coffee -reg -reg_method dynamic_msa -pool \
         -seq ${seqs} \
         -reg_tree ${guide_tree} \
         -reg_nseq ${bucket_size} \
         -dynamic ${dynamicX} \
         -reg_homoplasy \
         -dynamic_config ${dynamicConfig} \
         -outfile ${id}.dynamic_${params.db}_${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.with.${tree_method}.tree.aln
