#!/bin/bash

links_file=$1
output_dir=$2
file_end_pattern=$3

while read link; do
    if [[ $link =~ 'ftp://ftp.ncbi.nlm.nih.gov/geo/' ]]; then
	wget -nc -P $output_dir "$link/suppl/*$file_end_pattern"
    else
    	wget -nc -r $file_end_pattern -P $output_dir "$link"
    fi
done < $links_file

