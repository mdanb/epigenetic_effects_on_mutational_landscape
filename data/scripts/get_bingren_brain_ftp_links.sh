#!/bin/bash


esearch -db gds -query "GSE244618"  | efetch -format docsum | xtract -pattern DocumentSummary -if "title" -contains "ATAC-Seq" -element FTPLink > ../bingren_brain_ftp_links.txt
