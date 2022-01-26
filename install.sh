#!/usr/bin/bash
dockerURI="pdxfinder/pdx-liftover:vep_release98.3"
#printf "Pulling docker image at docker://%s in %s" "$dockerURI" "$(pwd)"
#singularity pull docker://"$dockerURI"
#singularity build "$dockerURI" 
singularity build --sandbox pdx-liftover-vep/ docker://"$dockerURI"

#printf "Installing database. Download is approximately 14G\n"
mkdir -p ./vepDB
cd ./vepDB || exit 1
mkdir -p homo_sapiens
cd ./homo_sapiens|| exit 1
curl -O ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
cd ..
curl -O ftp://ftp.ensembl.org/pub/release-98/variation/indexed_vep_cache/homo_sapiens_merged_vep_98_GRCh38.tar.gz
tar -zxcf homo_sapiens_merged_vep_98_GRCh38.tar.gz
