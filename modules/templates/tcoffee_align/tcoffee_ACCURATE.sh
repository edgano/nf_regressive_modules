t_coffee -seq $seqs -mode accurate -blast LOCAL -pdb_db ${params.database_path} $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_accurate.aln
