

run_pipeline() {
    local cancer_type="$1"
    echo $cancer_type
    echo "Running 2_Sorting_MutationFileSex_CancerType.py"
    python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types="$cancer_type"
    echo "Running 3_Intersect_paz_cancertypes.py"
    python3 3_Intersect_paz_cancertypes.py --cancer_types="$cancer_type" --increment_by=5 --max_samples=35 --subsample
    echo "Running 4_AssembleCout_paz_Cancergroup.py"
    python3 4_AssembleCout_paz_Cancergroup.py --cancer_types="$cancer_type" --subsampled
    echo "Running align_mutations_to_ranges.R"
    Rscript align_mutations_to_ranges.R --cancer_type="$cancer_type" --subsampled
}

run_pipeline "Lung-SCC"
#run_pipeline "Liver-HCC"
#run_pipeline "CNS-GBM"
#run_pipeline "ccRCC"
#run_pipeline "Skin-Melanoma"
#run_pipeline "CNS-Medullo"
