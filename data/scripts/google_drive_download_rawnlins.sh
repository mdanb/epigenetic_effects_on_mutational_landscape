while read -r line; do
    IFS=',' read -r filename link <<< "$line"
    gdown --fuzzy -O "../../data/bed_files/rawlins/${filename}.tsv.gz" "$link"
done < ../../data/bed_files/rawlins/links_and_names.txt
