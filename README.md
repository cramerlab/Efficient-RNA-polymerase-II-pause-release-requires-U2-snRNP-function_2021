# Efficient RNA-polymerase II pause release requires U2 snRNP function (2020)

### Fastq files and annotation used on this analysis are available from GEO series GSE148433.

#### Pre-analysis:   
TT-seq and RNA-seq reads aligned using STAR and filtered with Samtools. HTSeq  was used to calculate the read counts for different featues (please see https://github.com/cramerlab/TT-seq_analysis)


#### Annotation:
We included in our analysis only major isoforms with 70% or higher prevalenve per gene in both DMSO and Pla-B.
Major isoform annotation used on this analysis can be found in https://github.com/cramerlab/Efficient-RNA-polymerase-II-pause-release-requires-U2-snRNP-function_2020/Annotation. 
In order to create the major isoforms annotation:

#Download 

~/bash/salmon.sh

~/R_scripts/create_anno.R

~/R_scripts/major_isoform_selection.R


Intronless genes annotation used on this analysis can be found in https://github.com/cramerlab/Efficient-RNA-polymerase-II-pause-release-requires-U2-snRNP-function_2020/Annotation. 
In order to create the intronless genes annotatio:

~/Rscripts/intronless_annotation.R

#### Splicing ratio calculation:
1) Create Exon based spliced junction read count lists using bam files:

~/bash/exon_based_spliced_junction.R

2) Calculate the splicing ratio and plot the results:

~/R_scripts/splicing_ratio.R

