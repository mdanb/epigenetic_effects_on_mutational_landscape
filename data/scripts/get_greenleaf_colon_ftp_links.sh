esearch -db gds -query "GSE201336"  | efetch -format docsum | xtract -pattern DocumentSummary -element FTPLink | tail -n+4 > ../greenleaf_colon_ftp_links.txt
