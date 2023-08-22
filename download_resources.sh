#!/bin/bash

##................................................................###
## To set some colors ##
export COLOR_BLUE='\e[0;34m'
export COLOR_RED='\e[0;31m'
export COLOR_WHITE='\e[1;37m'

##Download URL

echo -e "${COLOR_BLUE} Do you want to Download hg38 from GENCODE.?"
read  -p "Enter [y/n]:" var1
echo -e "${COLOR_BLUE} Downloading Basic Annotation and FASTA file (Release 43). If not satisfied edit the link in the script to use latest release"

if [ $var1 == "y" ]; then
	mkdir resources;
	wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz \
  -P resources;
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.basic.annotation.gtf.gz \
  -P resources;
    gunzip resources/gencode.v43.primary_assembly.basic.annotation.gtf.gz;
	gunzip resources/GRCh38.primary_assembly.genome.fa.gz
	else
	echo -e "${COLOR_RED} GENCODE file download ABORTED";
fi

echo -e "${COLOR_WHITE} The resources will be downloaded in resources/ folder. \n"

wget  https://github.com/sortmerna/sortmerna/blob/master/data/rRNA_databases/silva-euk-18s-id95.fasta?raw=true -P resources
wget https://github.com/sortmerna/sortmerna/blob/master/data/rRNA_databases/silva-euk-28s-id98.fasta?raw=true -P resources
mv resources/silva-euk-18s-id95.fasta\?raw\=true resources/silva-euk-18s-id95.fasta
mv resources/silva-euk-28s-id98.fasta\?raw\=true resources/silva-euk-28s-id98.fasta

echo -e "${COLOR_WHITE} Silva database download completed..."