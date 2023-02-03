#!/bin/bash


esearch -db gds -query "GSE184462"  | efetch -format docsum | xtract -pattern DocumentSummary -if "title" -contains "ATAC-Seq" -element FTPLink > ../bingren_ftp_links.txt
