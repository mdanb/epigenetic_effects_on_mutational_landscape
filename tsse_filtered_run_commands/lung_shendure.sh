#!/bin/bash

Rscript ../create_tsse_filtered_count_overlaps.R --top_tsse_fragment_count_range="1000,10000,50000,100000,150000,250000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000,4000000,5000000,6000000,8000000,10000000,15000000,20000000,30000000,40000000,50000000" --dataset="shendure" --cell_types="Bronchiolar and alveolar epithelial cells" --files_pattern="lung" --cores=1
