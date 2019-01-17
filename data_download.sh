#!/bin/bash
# download and process databse file 
echo "download reference file"
mkdir database && cd database
cd database
wget http://141.2.194.197/rnaeditor_annotations/GRCH38.tar.gz
tar xvzf GRCH38.tar.gz && rm -f GRCH38.tar.gz
cd GRCH38 && mv ESP_filtered ESP.vcf
mkdir STAR_index
../../software/STAR --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --runThreadN 8
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
wget ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz && gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
cd ..
cd ..






