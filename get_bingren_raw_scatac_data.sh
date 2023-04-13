ftp_links=$(esearch -db gds -query "GSE184462" | 
			efetch -format docsum | 
			xtract -pattern DocumentSummary -if "title" -contains "ATAC-Seq" -element FTPLink)

for link in $ftp_links; do
	echo $link >> bingren_ftp_links.txt
done
