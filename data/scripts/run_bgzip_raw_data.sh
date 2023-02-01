#!/bin/bash

# argument = list of files
cat $1 | xargs -P8 -n1 ./bgzip_raw_data.sh
