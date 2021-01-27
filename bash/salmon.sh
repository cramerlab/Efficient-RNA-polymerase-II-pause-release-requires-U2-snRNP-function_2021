#!/bin/bash

#Download Refseq transcriptome
curl ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz -o Annotation/human_RefSeq_transc.fa.gz
gunzip Annotation/human_RefSeq_transc.fa.gz

#Defining paths
salmon_path=~/miniconda3/envs/salmon/bin/salmon
fastq_dir=~/Project/FastqFiles/

#Create index
${salmon_path} index -t /Annotation/human_RefSeq_transc.fa -i Annotation/human_RefSeq_transc_index

#Run salmon
#Fastq files example: control_rep1_read1.fastq, control_rep1_read2.fastq
mkdir ~/Salmon_quant
for f in $(ls ${fastq_dir}/*.fastq | rev | cut -c 12- | rev | uniq);
do file=`basename ${f}`
  ${salmon_path} quant -i Annotation/human_RefSeq_transc_index -l A -1 ${file}read1.fastq -2 ${file}read2.fastq -p 8 --validateMappings -o ~/Salmon_quant/${file}_quant
