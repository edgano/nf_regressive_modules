https://nextflow-io.github.io/patterns/index.html#_optional_input

t_coffee -reg -reg_method 3dcoffee_msa \
         -seq ${seqs} \
         -template_file ${template}
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_${bucket_size}.${align_method}.with.NO_TREE.tree.aln
