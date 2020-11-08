#!/bin/bash

## 
## add sufix
##
for file in *; do
    mv "$file" "${file}.fa"
done