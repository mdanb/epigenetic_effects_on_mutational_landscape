#!/bin/bash

Rscript ../create_tsse_filtered_count_overlaps.R --top_tsse_fragment_count_range="1000,10000,50000,100000,150000,250000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000,4000000,5000000,6000000" --dataset="shendure" --cell_types="" --files_pattern="brain" --cores=1
