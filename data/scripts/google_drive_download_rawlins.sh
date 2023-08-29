while read -r line; do
    IFS=',' read -r filename link <<< "$line"
    gdown --fuzzy -O "../bed_files/rawlins_fetal_lung_scATAC/${filename}.tsv.gz" "$link"
done < ../bed_files/rawlins_fetal_lung_scATAC/links_and_names.txt
