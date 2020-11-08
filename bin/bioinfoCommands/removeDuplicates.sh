#!/bin/bash

# The following pipline is composed of the following actions:
# in every row that starts with '>' put '@' at the end of the row
# replace '>' with '#'
# remove the break lines, replace '#' with a breakline, replace '@' with a tab
# sort the file according to the second column (sequences). the -u option keeps only one copy of each unique string.
# add '>' at the begining of each line
# sustitute tab (\t)  with a breakline (\n)
# remove the first line (it's a blank line) and save the result into $output

# To get only IDs of duplicated sequences with the same ID (Assuming duplicate records have identical IDs)
#   grep ">" dupl.fa | sort | uniq -d

FILES=./*.fa
for f in $FILES
do
    numDup="$(grep ">" $f | sort | uniq -d | wc -l)"
    #echo $numDup
    if [ "$numDup" -ne "0" ]; then                        # duplicates found
        echo "$numDup duplicates at $f";
        filename=${f%.[^.]*}
        sed '/^>/s/$/@/' < $f |
        sed 's/^>/#/' |\
        tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
        sort -u -f -k 2,2  |\
        sed -e 's/^/>/' |\
        tr '\t' '\n/' |\
        tail -n +2 > $filename'_noDup.fa'
        rm $filename'.fa'
    else
        echo "No duplicates at $f";
    fi
done