#!/bin/bash

cancer_types_file=$1
figure_dir="figures/models/XGB"
experiment_setting=$2

while read cancer_type; do
	from_coo $figure_dir/$cancer_type/$experiment_setting $figure_dir/$cancer_type/
done <$cancer_types_file
