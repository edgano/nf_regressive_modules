export blast_server_4_CLTCOFFEE=LOCAL
export protein_db_4_CLTCOFFEE=${params.database_path}
export VERBOSE_4_DYNAMIC=1

t_coffee -method psicoffee_msa \
         -tree ${guide_tree} \
         -seq ${seqs} \
         -cache ${params.blastOutdir} \
         -blast=LOCAL -protein_db=${params.database_path} \
         -outfile ${id}.prog.${align_method}.with.${tree_method}.tree.aln