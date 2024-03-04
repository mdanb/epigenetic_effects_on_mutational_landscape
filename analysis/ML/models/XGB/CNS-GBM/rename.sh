DIRECTORY="./"

# Loop through files matching the initial pattern
for file in ${DIRECTORY}scATAC_source_Bingren_Bingren_adult_brain_Greenleaf_brain_Shendure_cell_number_filter_100_tissues_to_consider_*; do
  if [[ $file =~ scATAC_source_Bingren_Bingren_adult_brain_Greenleaf_brain_Shendure_cell_number_filter_100_tissues_to_consider_adult_brain_brain_frontal_cortex_cerebrum_cerebellum_(.+) ]]; then
    # Capture the variable part of the filename after the pattern
    suffix="${BASH_REMATCH[1]}"
    
    # Construct the new filename with the correct order of tissues
    new_file="scATAC_source_Bingren_Bingren_adult_brain_Greenleaf_brain_Shendure_cell_number_filter_100_tissues_to_consider_adult_brain_frontal_cortex_cerebrum_brain_cerebellum_${suffix}"
    
    # Rename the file
    #mv "$file" "$new_file"
    echo "Renamed $file to $new_file"
  fi
done
