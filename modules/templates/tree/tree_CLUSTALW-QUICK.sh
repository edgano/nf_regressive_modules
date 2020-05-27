clustalw -infile=${seqs} -quicktree

mv ${id}.dnd  ${id}.${tree_method}.dnd
