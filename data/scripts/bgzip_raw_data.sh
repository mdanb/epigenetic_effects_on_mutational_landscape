#!/bin/bash

name=$(echo $1 | sed 's/...$//')
zcat $1 | sort -T /broad/hptmp/bgiotti -k1,1 -k2,2n | bgzip > ${name}.bgz
