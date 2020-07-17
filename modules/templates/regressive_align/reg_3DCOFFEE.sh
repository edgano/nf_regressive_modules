#! /bin/bash 

## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## 
## -- Structure to use optional INPUT files with templates -- ##
## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## 

declare template_filter=""
declare libs_filter=""

if [ $template == "input.4" ]; then
    \$template_filter=""
else
    \$template_filter=" -template_file $template"
fi

if [ $library == "input.5" ]; then
    \$libs_filter=""
else
    \$libs_filter=" -lib $library"
fi

## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##              https://nextflow-io.github.io/patterns/index.html#_optional_input
##         \$template_filter \$libs_filter \

t_coffee -reg -reg_method 3dcoffee_msa \
        -seq ${seqs} \
        -reg_nseq ${bucket_size} \
        -reg_homoplasy \
        \$template_filter \
        -outfile ${id}.reg_${bucket_size}.${align_method}.with.NO_TREE.tree.aln