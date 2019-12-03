#!/usr/bin/bash

echo "removing all vcf and ANN files"
find /nfs/nobackup/spot/mouseinformatics/pdx/data -name "*.tsv.vcf*" -exec rm {} \;
