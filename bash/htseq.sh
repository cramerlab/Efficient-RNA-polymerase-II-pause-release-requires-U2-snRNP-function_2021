samtools_path=~/samtools/bin/samtools
genome_anno_gtf=~/Human/Annotation/refseq_anno.gtf
spikeins_anno_gtf=~/Human/Annotation/spikeins.gtf
bam_dir=~/Project/BamFiles
counts_dir=~/Project/CountFiles
mkdir ${counts_dir}
 
#${star_path}  sort -n ${bam_dir}/.${bam} | ${star_path} view -h - | gzip > \${intermediatefolder}$2.namesorted.sam.gz

mkdir ${counts_dir}/Refseq
python -m HTSeq.scripts.count -m intersection-strict -t transcript -i gene_id --stranded=yes ${bam_dir}/.${bam} \$anno > ${counts_dir}/Refseq/.${tabular}
mkdir ${counts_dir}/RefseqAntissense
python -m HTSeq.scripts.count -m intersection-strict -t transcript -i gene_id --stranded=reverse ${bam_dir}/.${bam} \$anno > ${counts_dir}/RefseqAntissense/.${tabular}

mkdir ${counts_dir}/SpikeIns
python -m HTSeq.scripts.count -m intersection-strict -t transcript -i gene_id --stranded=yes ${bam_dir}/.${bam} \$anno > ${counts_dir}/SpikeIns/.${tabular}
mkdir ${counts_dir}/SpikeInsAntissense
python -m HTSeq.scripts.count -m intersection-strict -t transcript -i gene_id --stranded=reverse ${bam_dir}/.${bam} \$anno > ${counts_dir}/SpikeInsAntissense/.${tabular}
