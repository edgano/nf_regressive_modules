export blast_server_4_CLTCOFFEE=LOCAL
export protein_db_4_CLTCOFFEE=${params.database_path}
export VERBOSE_4_DYNAMIC=1

t_coffee -method psicoffee_msa \
         -tree ${guide_tree} \
         -seq ${seqs} \
         -cache ${params.blastOutdir} \
         -psitrim 100 -psiJ 3  -prot_min_cov 90  -prot_max_sim 100         -prot_min_sim 0 \
         -blast=LOCAL -protein_db=${params.database_path} \
         -outfile ${id}.prog.${align_method}.with.${tree_method}.tree.aln