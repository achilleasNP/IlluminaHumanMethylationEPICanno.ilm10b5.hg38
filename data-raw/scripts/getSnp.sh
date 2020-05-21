#!/bin/bash -e
mkdir  -p files  || exit 1
cd files || exit 1
for i in {141,142,144,146,147,150,151}; do
    filename="snp${i}Common.txt.gz"
    if [[ ! -e ${filename} ]]; then
        curl -O "ftp://hgdownload.soe.ucsc.edu/apache/htdocs/goldenPath/hg38/database/${filename}"
    fi
    gunzip -cd "snp${i}Common.txt.gz" | cut -f2,3,4,5,7,8,12,25 | gzip -c > "snp${i}Common_small.txt.gz"
done

