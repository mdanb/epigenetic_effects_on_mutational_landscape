#!/bin/bash

Rscript ../create_tsse_filtered_count_overlaps.R --top_tsse_fragment_count_range="1000,10000,50000,100000,150000,250000,300000,322928,322929,400000,500000,600000" --dataset="bing_ren" --cell_types="Foveolar Cell" --files_pattern="esophagus_muscularis" --cores=1
