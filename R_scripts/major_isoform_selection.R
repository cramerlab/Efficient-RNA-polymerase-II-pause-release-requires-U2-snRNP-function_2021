#libraries needed
library(biomaRt)
library(dplyr)
library(
library(rtracklayer)
library(GenomicRanges)

 #notin function 
`%notin%` <- Negate(`%in%`)
  
  
#COMPARE DMSO1h vs PlaB1h
DMSO1h_1=read.table("../salmon/DMSO1h_1_quant.sf", header = T)
DMSO1h_2=read.table("../salmon/DMSO1h_2_quant.sf", header = T)
PlaB1h_1=read.table("../salmon/PlaB1h_1_quant.sf", header = T)
PlaB1h_2=read.table("../salmon/PlaB1h_2_quant.sf", header = T)


Quant=cbind(DMSO1h_1[,4],DMSO1h_2[,4],PlaB1h_1[,4], PlaB1h_2[,4])
dim(Quant)

colnames(Quant)=c("DMSO1h_1", "DMSO1h_2", "PlaB1h_1", "PlaB1h_2")
rownames(Quant)=PlaB1h_1$Name

#Remove rows with all values inferior to 0.5 TPM 
Quant=Quant[apply(Quant[,-1], 1, function(x) !all(x<0.5)),]
dim(Quant)

#Remove model RefSeqs XR, XM
Quant=Quant[- grep("XR", rownames(Quant)),]
Quant=Quant[- grep("XM", rownames(Quant)),]
dim(Quant) 

#Add column with mean 
Quant=cbind(Quant, "mean_PlaB"= rowMeans(Quant[,c("PlaB1h_1", "PlaB1h_2")]))
Quant=cbind(Quant, "mean_DMSO"= rowMeans(Quant[,c("DMSO1h_1", "DMSO1h_2")]))
Quant=cbind(Quant, "mean_all"= rowMeans(Quant))

#Load refseq annotation file
load("human.refseq.anno.April19.RData")
                  
#Remove the version number, because I dont have it on the biomart anno
rownames(Quant)=sub('\\..*', '', rownames(Quant))
human.refseq.anno$transcript_id=sub('\\..*', '', human.refseq.anno$transcript_id)

                 
                  
load("ensembl.biomart.refseq.2.ensembl.gene.id.mapping.RData")


human.refseq.anno[,"gene_name"]=ensembl.biomart.refseq.2.ensembl.gene.id.mapping[human.refseq.anno[,"transcript_id"]]

                  
refseq_salmon=human.refseq.anno[,c("gene_name", "transcript_id")]
transcripts = rownames(Quant) 
refseq_salmon=refseq_salmon[which(refseq_salmon$transcript_id %in% transcripts),]
#length(unique(refseq_salmon$gene_id)) 
#there was two diffferent gene_ids for the same gene. Thats why now I am using gene names

length(unique(refseq_salmon$gene_name))
length(unique(refseq_salmon$transcript_id)) 
library(dplyr)
refseq_salmon=distinct(refseq_salmon)
length(unique(refseq_salmon$gene_name))
length(unique(refseq_salmon$transcript_id)) 
length(setdiff(rownames(Quant), refseq_salmon$transcript_id )) 

dif=setdiff(rownames(Quant), refseq_salmon$transcript_id )
length(dif) 

Quant <- Quant[ ! rownames(Quant) %in% dif, ]
dim(Quant)
refseq_salmon=cbind(refseq_salmon, Quant[as.vector(refseq_salmon$transcript_id),])
dim(refseq_salmon)
#refseq_salmon=refseq_salmon[apply(refseq_salmon[,3:6], 1, function(x) !all(x==0)),]

refseq_salmon[refseq_salmon<0.001] <- NA
refseq_salmon<-refseq_salmon[complete.cases(refseq_salmon),]
length(unique(refseq_salmon$gene_name))

length(unique(refseq_salmon$transcript_id)) 
    
load("human.refseq.anno.April19.RData")

# Calculate the percentage of each isoform                  
refseq_salmon_filtered=refseq_salmon %>%
  group_by(gene_name) %>%
  mutate(percDMSO = mean_DMSO/sum(mean_DMSO))  %>%
  mutate(percPlaB = mean_PlaB/sum(mean_PlaB))  %>%
  mutate(percall = mean_all/sum(mean_all))  
                  
##Select isoforms with at least 70% prevalance and with mean TMP > 0.5 on DMSO (control)
refseq_salmon_filtered=refseq_salmon_filtered %>%
  group_by(gene_name) %>%
  filter(mean_all == max(mean_all))

#remove isoforms with less than 70% in DMSO and PlaB
refseq_salmon_filtered=refseq_salmon_filtered[which(refseq_salmon_filtered[,"percDMSO"] >= 0.7),] 
refseq_salmon_filtered=refseq_salmon_filtered[which(refseq_salmon_filtered[,"percPlaB"] >= 0.7),] 

#remove isoforms with less than 0.5 mean TPM on DMSO samples
refseq_salmon_filtered=refseq_salmon_filtered[which(refseq_salmon_filtered[,"mean_DMSO"] >= 0.5),]

length(unique(refseq_salmon_filtered$gene_name)) 
length(unique(refseq_salmon_filtered$transcript_id)) 

isoforms=refseq_salmon_filtered$transcript_id
                  
#Generate the annotation file using refseq
load("human.refseq.extended.April19.RData")
human.refseq.anno.TR = human.refseq.extended[which(human.refseq.extended[,"type"] == "transcript"),]
human.refseq.anno.TR$transcript_id=sub('\\..*', '', human.refseq.anno.TR$transcript_id)
human.refseq.major.TR = human.refseq.anno.TR[which(human.refseq.anno.TR[,"transcript_id"] %in% isoforms),]
human.refseq.major.TR=human.refseq.major.TR[,c("chr", "strand", "type", "start", "end", "gene_name", "transcript_id", "length", "id"),]
length(human.refseq.major.TR$transcript_id)
length(unique(human.refseq.major.TR$transcript_id)) 
                  

#Remove chr X, Y and M
human.refseq.major.TR=human.refseq.major.TR[which(human.refseq.major.TR$chr %notin% c("chrX","chrY","chrM")),]
rm(human.refseq.major.TR)
length(unique(human.refseq.major.TR$transcript_id))
#save(human.refseq.major.TR_Oct2020, file= file.path("AnnotationObjects","human.refseq.major.TR_Oct2020.RData"))

#Find and remove overlapping major isoforms                  
gr = makeGRangesFromDataFrame(human.refseq.major.TR, keep.extra.columns = TRUE)
hits <- findOverlaps(gr, type = "any", select = "all", ignore.strand=FALSE, drop.self=TRUE, drop.redundant=FALSE)
ovpairs <- Pairs(gr, gr, hits=hits)
pint <- pintersect(ovpairs, ignore.strand=TRUE)
ids_with_overlaps=as.character(names(pint))

human.refseq.major.TR_overlapped=human.refseq.major.TR[which(human.refseq.major.TR$id %in% ids_with_overlaps),]
TRs_overlapped=human.refseq.major.TR_overlapped$transcript_id
length(TRs_overlapped)
human.refseq.major.TR=human.refseq.major.TR[which(human.refseq.major.TR$transcript_id %notin% TRs_overlapped),]
length(human.refseq.major.TR$transcript_id)
length(unique(human.refseq.major.TR$transcript_id))
                  
###Find and remove overlapping annotated transcripts
#Extended with major isoforms selected
human.refseq.extended$transcript_id=sub('\\..*', '', human.refseq.extended$transcript_id)
human.refseq.extended.tr=human.refseq.extended[which(human.refseq.extended$type == "transcript"),]
human.refseq.extended.tr.major=human.refseq.extended.tr[which(human.refseq.extended.tr$transcript_id %in% human.refseq.major.TR$transcript_id),]

#Extended without the genes from the major isoforms
human.refseq.extended.tr.not.major=human.refseq.extended.tr[which(human.refseq.extended.tr$gene_name %notin% human.refseq.extended.tr.major$gene_name),]

#merge extended with major isoforms + extended without genes for major isoforms
human.refseq.extended.merge=rbind(human.refseq.extended.tr.major,human.refseq.extended.tr.not.major)

gr = makeGRangesFromDataFrame(human.refseq.extended.merge, keep.extra.columns = TRUE)
hits <- findOverlaps(gr, type = "any", minoverlap=0, select = "all", ignore.strand=FALSE, drop.self=TRUE, drop.redundant=FALSE)

ovpairs <- Pairs(gr, gr, hits=hits)
pint <- pintersect(ovpairs, ignore.strand=FALSE)

ids_with_overlaps=pint$transcript_id

human.refseq.major.TR=human.refseq.major.TR[which(human.refseq.major.TR$transcript_id %notin% ids_with_overlaps),]
length(human.refseq.major.TR$transcript_id)
                  
human.refseq.major.TR_overlap=human.refseq.major.TR[which(human.refseq.major.TR$transcript_id %in% ids_with_overlaps),]
length(human.refseq.major.TR_overlap$transcript_id)


save(human.refseq.major.TR,  file= "Annotation/human.refseq.major.TR.RData")
export(human.refseq.major.TR,con = "Annotation/human.refseq.major.TR.gtf")
                  
                  
