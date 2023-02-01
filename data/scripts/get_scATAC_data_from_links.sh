#!/bin/bash

links_file=$1
output_dir=$2

while read link; do
    if [[ $link =~ 'ftp://ftp.ncbi.nlm.nih.gov/geo/' ]]; then
	file_end_pattern=$3
	wget -nc -P $output_dir "$link/suppl/*$file_end_pattern"
    else
    	wget -nc -P $output_dir $link
    fi
done < $links_file

