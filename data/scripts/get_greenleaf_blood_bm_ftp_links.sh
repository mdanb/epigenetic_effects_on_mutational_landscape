esearch -db gds -query "GSE139369" | efetch -format docsum | xtract -pattern DocumentSummary -if "title" -contains "ATAC" -element FTPLink > ../greenleaf_blood_bm_ftp_links.txt
