#!/bin/bash

Rscript ../create_tsse_filtered_count_overlaps.R --top_tsse_fragment_count_range="1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000" --dataset="tsankov" --cell_types="AT1;AT2;B_cells;Ciliated;Endothelial;Fibroblasts;Immune;Mesothelium;Secretory;SmoothMuscle" --files_pattern="RPL" --cores=4
