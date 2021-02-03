load("ProcessedData/exon.based.spliced.junction.read.count.list_all.non.single.RData")

#sequencing.depth calculation can be found in https://github.com/cramerlab/TT-seq_analysis/
load("ProcessedData/sequencing.depth.RData") 

control_1 <- exon.based.spliced.junction.read.count.list$control1.bam[,c(6,12:19)] 
control_2 <- exon.based.spliced.junction.read.count.list$control2.bam[,c(6,12:19)]
treatment_1 <- exon.based.spliced.junction.read.count.list$treatment1.bam[,c(6,12:19)] 
treatment_2 <- exon.based.spliced.junction.read.count.list$treatment2.bam[,c(6,12:19)]

treatment_1 = treatment_1[which(treatment_1$five.prime.total>=30),]
dim(treatment_1) 
treatment_2=treatment_2[which(rownames(treatment_2) %in% rownames(treatment_1)),]
dim(AMOU2_2)
control_1=CTR_1[which(rownames(control_1) %in% rownames(treatment_1)),]
control_2=CTR_2[which(rownames(control_2) %in% rownames(treatment_1)),]
dim(control_2)

#5'SS
control_5prime = cbind(control_1[c("five.prime.spliced", "five.prime.unspliced")], control_2[c("five.prime.spliced","five.prime.unspliced")])
colnames(control_5prime) = c("control_1_spliced","control_1_unspliced", "control_2_spliced", "control_2_unspliced")

treatment_5prime = cbind(treatment_1[c("five.prime.spliced", "five.prime.unspliced")], treatment_2[c("five.prime.spliced","five.prime.unspliced")])
colnames(treatment_5prime) = c("treatment_1_spliced","treatment_1_unspliced", "treatment_2_spliced", "treatment_2_unspliced")

countData=cbind(control_5prime, treatment_5prime)

colData = data.frame("condition" =c("control","control","control","control","treatment","treatment","treatment","treatment"), 
                     "assay" = c("spliced","unspliced","spliced","unspliced","spliced","unspliced","spliced","unspliced") ,row.names = colnames(countData))
dds = DESeqDataSetFromMatrix(countData = countData, 
                             colData = colData,
                             design =  ~ assay + condition + assay:condition)

dds$condition = relevel(dds$condition,ref="control")

names(sequencing.depth)=c("treatment_1", "control_1","control_2","treatment_2")
normFactors = matrix(rep(sequencing.depth[substr(names(countData),1,7)],each=nrow(countData)),nrow=nrow(countData),ncol=ncol(countData))
normalizationFactors(dds) = normFactors/mean(normFactors)

dds = dds[rowSums(counts(dds)) > 1,]

##Affected
dds = DESeq(dds)
res_5SS <- results(dds, contrast=c("condition","treatment", "control"), altHypothesis = "less") 
summary(res_5SS)
save(res_5SS,file="res_5SS_spikein_exons.RData")
plotMA(res_5SS, main="Splicing 5'SS (Exons)", ylim=c(-3,3))

#Unnafected
ddsNoPrior <- DESeq(dds, betaPrior=FALSE)
  
unchange_5SS=results(ddsNoPrior, contrast=c("condition","treatment", "control"), lfcThreshold=1.5, altHypothesis = "lessAbs")
summary(unchange_5SS)

save(unchange_5SS,file="unchange_5SS_spikein_exons.RData")

plotMA(unchange_5SS, main="Unchanged 5'SS (Exons), altHypothesis=leAMOU2bs",ylim=c(-3,3))

#3'SS
control_3prime = cbind(control_1[c("three.prime.spliced", "three.prime.unspliced")], control_2[c("three.prime.spliced","three.prime.unspliced")])
colnames(control_3prime) = c("control_1_spliced","control_1_unspliced", "control_2_spliced", "control_2_unspliced")

treatment_3prime = cbind(treatment_1[c("three.prime.spliced", "three.prime.unspliced")], treatment_2[c("three.prime.spliced","three.prime.unspliced")])
colnames(treatment_3prime) = c("treatment_1_spliced","treatment_1_unspliced", "treatment_2_spliced", "treatment_2_unspliced")

countData=cbind(CTR_3prime, AMOU2_3prime)

colData = data.frame("condition" =c("control","control","control","control","treatment","treatment","treatment","treatment"), 
                     "assay" = c("spliced","unspliced","spliced","unspliced","spliced","unspliced","spliced","unspliced") ,row.names = colnames(countData))
dds = DESeqDataSetFromMatrix(countData = countData, 
                             colData = colData,
                             design =  ~ assay + condition + assay:condition)

dds$condition = relevel(dds$condition,ref="control")

normalizationFactors(dds) = normFactors/mean(normFactors)

dds = dds[rowSums(counts(dds)) > 1,]

##Affected
dds = DESeq(dds)
#res = results(dds)
res_3SS <- results(dds, contrast=c("condition","treatment", "control"), altHypothesis = "less") 
summary(res_3SS)
save(res_3SS,file="res_3SS_spikein_exons.RData")
plotMA(res_5SS, main="Splicing 3'SS (Exons)", ylim=c(-3,3))

#Unnafected
ddsNoPrior <- DESeq(dds, betaPrior=FALSE)
  
unchange_3SS=results(ddsNoPrior, contrast=c("condition","treatment", "control"), lfcThreshold=1.5, altHypothesis = "lessAbs")
summary(unchange_3SS)

save(unchange_3SS,file="unchange_3SS_spikein_exons.RData")

plotMA(unchange_3SS, main="Unchanged 3'SS (Exons), altHypothesis=less",ylim=c(-3,3))

#Select affected exons
##Filter for padj < 0.05 & log2FoldChange < -0.5
five.prime.exons= subset(res_5SS,padj < 0.05 & log2FoldChange < -1.5)
summary(five.prime.exons)
five.prime.exons= rownames(five.prime.exons)

##Filter for padj < 0.05
three.prime.exons=subset(res_3SS,padj < 0.05 & log2FoldChange < -1.5)
summary(three.prime.exons)
three.prime.exons=rownames(three.prime.exons)

##Transcripts containing the affected 1st exon 5'SS and/or 2nd exon 3'SS
load("Annotation/human.refseq.major.isoform.exon.RData")
first.exon=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$exon_number==1),]
TR_five.prime.first=first.exon[which(first.exon$id %in% five.prime.exons),]
TR_five.prime.first.exon=unique(TR_five.prime.first$id)
TR_five.prime.first=unique(TR_five.prime.first$transcript_id)
length(TR_five.prime.first) 

sec.exon=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$exon_number==2),]
TR_three.prime.sec=sec.exon[which(sec.exon$id %in% three.prime.exons),]
TR_three.prime.sec.exon=unique(TR_three.prime.sec$id)
TR_three.prime.sec=unique(TR_three.prime.sec$transcript_id)
length(TR_three.prime.sec)

TR.exon=union(TR_five.prime.first, TR_three.prime.sec)
length(TR.exon) 

save(TR.exon, file= "TR.exon.RData")

#Select unaffected exons
five.prime.exons.un=subset(unchange_5SS,padj < 0.05 & log2FoldChange > -1.5 & log2FoldChange < 1.5)
five.prime.exons.un=rownames(five.prime.exons.un)
three.prime.exons.un=subset(unchange_3SS,padj < 0.05 & log2FoldChange > -1.5 & log2FoldChange < 1.5)
three.prime.exons.un=rownames(three.prime.exons.un)

##Transcripts containing the unaffected 1st exon 5'SS and/or 2nd exon 3'SS
first.exon=human.refseq.major.isoform.exon_Oct2020[which(human.refseq.major.isoform.exon$exon_number==1),]
TR_five.prime.un.first=first.exon[which(first.exon$id %in% five.prime.exons.un),]
TR_five.prime.un.first.exon=unique(TR_five.prime.un.first$id)
TR_five.prime.un.first=unique(TR_five.prime.un.first$transcript_id)
length(TR_five.prime.un.first) 

sec.exon=human.refseq.major.isoform.exon_Oct2020[which(human.refseq.major.isoform.exon$exon_number==2),]
TR_three.prime.un.sec=sec.exon[which(sec.exon$id %in% three.prime.exons.un),]
TR_three.prime.un.sec.exon=unique(TR_three.prime.un.sec$id)
TR_three.prime.un.sec=unique(TR_three.prime.un.sec$transcript_id)
length(TR_three.prime.un.sec)

TR.un.exon=union(TR_five.prime.un.first, TR_three.prime.un.sec)
length(TR.un.exon) 
length(intersect(TR.exon, TR.un.exon))

save(TR.un.exon, file= "TR.un.exon.RData")
