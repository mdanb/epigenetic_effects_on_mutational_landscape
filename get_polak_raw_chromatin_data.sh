#!/bin/bash

ftp_links=$(esearch -db gds -query "GSE18927"  | 
	    efetch -format docsum | 
	    xtract -pattern DocumentSummary -if "title" -contains "Chromatin accessibility" -element FTPLink)

for link in $ftp_links; do                                                      
    wget -nc -P raw_data "$link/suppl/*.bed.gz"                                                 
done
