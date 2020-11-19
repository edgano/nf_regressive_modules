#! /bin/bash 

declare compressFlag=" "

if $params.compressAZ ; then
    compressFlag=" -output fastaz_aln"
fi

{ time -p t_coffee -reg -reg_method probcons_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         \$compressFlag \
         -outfile ${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln 2> tcoffee.stderr ; } 2> time.txt
