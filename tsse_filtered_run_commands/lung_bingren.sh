#!/bin/bash

Rscript ../create_tsse_filtered_count_overlaps.R --top_tsse_fragment_count_range="1000,10000,50000,100000,150000,250000,300000,400000,500000,600000" --dataset="bing_ren" --cell_types="Alveolar Type 2 (AT2) Cell" --files_pattern="lung_SM"
