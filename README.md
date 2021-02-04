# Efficient RNA-polymerase II pause release requires U2 snRNP function (2021)

### Fastq files and annotation used on this analysis are available from GEO series GSE148433.

#### Pre-analysis:   
TT-seq and RNA-seq reads were aligned using STAR and filtered with Samtools. HTSeq  was used to calculate the read counts for different featues (please see https://github.com/cramerlab/TT-seq_analysis)

#### Annotation:
We included in our analysis only major isoforms with 70% or higher prevalenve per gene in both DMSO and Pla-B.
In order to create the major isoforms annotation:  

##### Create general annotation:
Download ncbi Refseq hg38 genome assembly from https://genome.ucsc.edu/cgi-bin/hgTables   
(clade: Mammal, genome: Human,  assembly: Dex.2013 GRCh38/h38, group: Gene and Gene Predictions, track: NCBI RefSeq, table: RefSeq Curated)   
R_scripts/create_anno.R   
R_scripts/create_extended.R 

##### Create major isoform annotation:
bash/salmon.sh  
R_scripts/major_isoform_selection.R. 

##### Create intronless genes annotation:
R_scripts/intronless_annotation.R   

#### Splicing ratio calculation:
1) Create Exon based spliced junction read count lists using bam files:  
R_scripts/exon_based_spliced_junction.R   
2) Calculate the splicing ratio and plot the results:  
R_scripts/splicing_ratio.R

#### Identification of splicing-affected and unaffected transcripts:
R_scripts/spicing_affected_unaffected.R

#### Calculating elongation velocity:
R_scripts/coverage_rle_tracks.R  
R_scripts/coverage.R  
R_scripts/RNA_amount_per_cell.R  
R_scripts/elong_velocity.R  
