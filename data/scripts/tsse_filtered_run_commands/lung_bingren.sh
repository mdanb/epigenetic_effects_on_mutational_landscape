#!/bin/bash

Rscript ../create_tsse_filtered_count_overlaps.R --top_tsse_fragment_count_range="1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,15000000,20000000" --dataset="bing_ren" --cell_types="Alveolar Type 2 (AT2) Cell,Club Cell" --files_pattern="lung_SM" --cores=1
