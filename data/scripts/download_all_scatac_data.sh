dirs=("bingren_scATAC" "bingren_adult_brain" "JShendure_scATAC" "Tsankov_scATAC" "greenleaf_brain_scATAC" "greenleaf_pbmc_bm_scATAC" "greenleaf_colon_scATAC" "rawlins_fetal_lung_scATAC")

for dir in "${dirs[@]}"; do
	mkdir -p ../bed_files/$dir
done

# Bingren Atlas
./get_scATAC_data_from_links.sh ../bingren_ftp_links.txt ../bed_files/bingren_scATAC *gz > /dev/null 2>&1 &

# Bingren Brain
./get_scATAC_data_from_links.sh ../bingren_brain_ftp_links.txt ../bed_files/bingren_adult_brain *gz > /dev/null 2>&1 & 

# Greenleaf Blood / Bone Marrow 
./get_scATAC_data_from_links.sh ../greenleaf_blood_bm_ftp_links.txt ../bed_files/greenleaf_pbmc_bm_scATAC *gz > /dev/null 2>&1 &

# Greenleaf Colon
./get_scATAC_data_from_links.sh ../greenleaf_colon_ftp_links.txt ../bed_files/greenleaf_colon_scATAC *gz > /dev/null 2>&1 &

# Greenleaf Fetal Brain
./get_scATAC_data_from_links.sh ../greenleaf_fetal_brain_links.txt ../bed_files/greenleaf_brain_scATAC *tsv > /dev/null 2>&1 & 

# Rawlins Fetal Lung
./google_drive_download_rawlins.sh > /dev/null 2>&1 &

# Shendure Fetal Atlas
./get_scATAC_data_from_links.sh ../shendure_link.txt ../bed_files/JShendure_scATAC *gz > /dev/null 2>&1 & 

