#!/bin/bash 

## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## 
## -- Structure to use optional INPUT files with templates -- ##
## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## 

declare template_filter=""
declare libs_filter=""

if [ $template == "input.3" ]; then
    \$template_filter=""
else
    \$template_filter="-template_file $template"
fi

if [ $library == "input.4" ]; then
    \$libs_filter=""
else
    \$libs_filter="-lib $library"
fi

## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## 

t_coffee -seq $seqs -mode mcoffee \$template_filter \$libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_mcoffee.aln
