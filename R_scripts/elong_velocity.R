load("AnnotationObjects/human.refseq.major.isoform.exon.RData")
load("AnnotationObjects/human.refseq.major.isoform.intron.RData")
load("AnnotationObjects/protein.coding.RData")

#sequencing.depth, cross.contamination and expressed.TR.5.50 calculation can be found in https://github.com/cramerlab/TT-seq_analysis/
load("ProcessedData/TTseq/sequencing.depth.RData")
load("ProcessedData/TTseq/cross.contamination.RData")
load("ProcessedData/TTseq/expressed.TR.5.50.RData")
load("ProcessedData/TTseq/protein.coding.RData")

#size.factors.cerevisiae calculation can be found in https://github.com/cramerlab/mNET-seq_analysis/
load("ProcessedData/mNETseq/size.factors.cerevisiae.RData")

load("ProcessedData/TTseq/conversion.factor.to.amount.per.cell.RData")
s=383 #calibration factor s was calibrated to resemble pause duration of previous experiments in K562 cells (Gressel et al., 2019; Gressel et al., 2017) .

human.refseq.major.isoform.exon=as.data.frame(human.refseq.major.isoform.exon)

human.refseq.major.isoform.exon.coverage.antisense.corrected.TTseq = get(load("ProcessedData/TTseq/human.refseq.major.isoform.exon.coverage.antisense.corrected.RData"))
human.refseq.major.isoform.intron.coverage.antisense.corrected.TTseq = get(load("ProcessedData/TTseq/human.refseq.major.isoform.intron.coverage.antisense.corrected.RData"))

human.refseq.major.isoform.exon.coverage.antisense.corrected.mNETseq = get(load("ProcessedData/mNETseq/human.refseq.major.isoform.exon.coverage.antisense.corrected.RData"))
human.refseq.major.isoform.intron.coverage.antisense.corrected.mNETseq = get(load("ProcessedData/mNETseq/human.refseq.major.isoform.intron.coverage.antisense.corrected.RData"))

#CONTROL
#TT-seq: control_1, control_2
#RNA-seq: RNAseq_control_1, RNAseq_control_2
#mNET-seq: control_1, control_2
# on first exons (non-single)
human.refseq.major.isoform.exon.first.non.single = human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon[,"exon_order"] == "first"),]
human.refseq.major.isoform.exon.first.non.single = human.refseq.major.isoform.exon.first.non.single[which(as.character(human.refseq.major.isoform.exon.first.non.single[,"transcript_id"]) %in% expressed.TR.5.50),]
human.refseq.major.isoform.exon.first.non.single = human.refseq.major.isoform.exon.first.non.single[which(as.character(human.refseq.major.isoform.exon.first.non.single[,"transcript_id"]) %in% protein.coding),]

exon.first.non.single = as.character(human.refseq.major.isoform.exon.first.non.single[,"id"])	
exon.first.non.single.list = tapply(as.character(human.refseq.major.isoform.exon.first.non.single[,"id"]),INDEX = as.character(human.refseq.major.isoform.exon.first.non.single[,"transcript_id"]),identity)

TTseq.human.refseq.major.isoform.exon.first.non.single.coverage.antisense.corrected = sapply(exon.first.non.single.list,function(x){sum(apply(t(t((t(t(human.refseq.major.isoform.exon.coverage.antisense.corrected.TTseq[x,c("control_1.bam","control_2.bam")])/sequencing.depth[c("control_1.tabular","control_2.tabular")])) - t(cross.contamination[c("control_1.tabular","control_2.tabular")]*t((t(t(human.refseq.major.isoform.exon.coverage.antisense.corrected.TTseq[x,c("RNAseq_control_1.bam","RNAseq_control_2.bam")])/sequencing.depth[c("RNAseq_control_1.tabular","RNAseq_control_2.tabular")])))))/(1 - cross.contamination[c("control_1.tabular","control_2.tabular")])),1,function(x){x[x < 0] = 0;sum(x)}))})

mNETseq.human.refseq.major.isoform.exon.first.non.single.coverage.antisense.corrected = sapply(exon.first.non.single.list,function(x){sum(human.refseq.major.isoform.exon.coverage.antisense.corrected.mNETseq[x,c("control_1.bam","control_2.bam")]/size.factors[c("control_1.tabular","control_2.tabular")])})

elongation.rate.exons.first.non.single = (TTseq.human.refseq.major.isoform.exon.first.non.single.coverage.antisense.corrected/mNETseq.human.refseq.major.isoform.exon.first.non.single.coverage.antisense.corrected)*(mNETseq.size.factor/(10*conversion.factor.to.amount.per.cell))

save(elongation.rate.exons.first.non.single,file="ProcessedData/elongation.rate.exons.first.control.RData")

# on first introns (including-single)
human.refseq.major.isoform.intron.first.non.single = human.refseq.major.isoform.intron[which(human.refseq.major.isoform.intron[,"intron_order"] == "first" | human.refseq.major.isoform.intron[,"intron_order"] == "1"),]
human.refseq.major.isoform.intron.first.non.single = human.refseq.major.isoform.intron.first.non.single[which(as.character(human.refseq.major.isoform.intron.first.non.single[,"transcript_id"]) %in% expressed.TR.5.50),]
human.refseq.major.isoform.intron.first.non.single = human.refseq.major.isoform.intron.first.non.single[which(as.character(human.refseq.major.isoform.intron.first.non.single[,"transcript_id"]) %in% protein.coding),]

intron.first.non.single = as.character(human.refseq.major.isoform.intron.first.non.single[,"id"])	
intron.first.non.single.list = tapply(as.character(human.refseq.major.isoform.intron.first.non.single[,"id"]),INDEX = as.character(human.refseq.major.isoform.intron.first.non.single[,"transcript_id"]),identity)

TTseq.human.refseq.major.isoform.intron.first.non.single.coverage.antisense.corrected = sapply(intron.first.non.single.list,function(x){sum(apply(t(t((t(t(human.refseq.major.isoform.intron.coverage.antisense.corrected.TTseq[x,c("control_1.bam","control_2.bam")])/sequencing.depth[c("control_1.tabular","control_2.tabular")])) - t(cross.contamination[c("control_1.tabular","control_2.tabular")]*t((t(t(human.refseq.major.isoform.intron.coverage.antisense.corrected.TTseq[x,c("RNAseq_control_1.bam","RNAseq_control_2.bam")])/sequencing.depth[c("RNAseq_control_1.tabular","RNAseq_control_2.tabular")])))))/(1 - cross.contamination[c("control_1.tabular","control_2.tabular")])),1,function(x){x[x < 0] = 0;sum(x)}))})

mNETseq.human.refseq.major.isoform.intron.first.non.single.coverage.antisense.corrected = sapply(intron.first.non.single.list,function(x){sum(human.refseq.major.isoform.intron.coverage.antisense.corrected.mNETseq[x,c("control_1.bam","control_2.bam")]/size.factors[c("control_1.tabular","control_2.tabular")])})

elongation.rate.introns.first.non.single = (TTseq.human.refseq.major.isoform.intron.first.non.single.coverage.antisense.corrected/mNETseq.human.refseq.major.isoform.intron.first.non.single.coverage.antisense.corrected)*(mNETseq.size.factor/(10*conversion.factor.to.amount.per.cell))

save(elongation.rate.introns.first.non.single,file="ProcessedData/elongation.rate.introns.first.control.RData")

#TREATMENT
#TT-seq: treatment_1, treatment_2
#RNA-seq: RNAseq_treatment_1, RNAseq_treatment_2
#mNET-seq: treatment_1, treatment_2
# on first exons (non-single)
human.refseq.major.isoform.exon.first.non.single = human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon[,"exon_order"] == "first"),]
human.refseq.major.isoform.exon.first.non.single = human.refseq.major.isoform.exon.first.non.single[which(as.character(human.refseq.major.isoform.exon.first.non.single[,"transcript_id"]) %in% expressed.TR.5.50),]
human.refseq.major.isoform.exon.first.non.single = human.refseq.major.isoform.exon.first.non.single[which(as.character(human.refseq.major.isoform.exon.first.non.single[,"transcript_id"]) %in% protein.coding),]

exon.first.non.single = as.character(human.refseq.major.isoform.exon.first.non.single[,"id"])	
exon.first.non.single.list = tapply(as.character(human.refseq.major.isoform.exon.first.non.single[,"id"]),INDEX = as.character(human.refseq.major.isoform.exon.first.non.single[,"transcript_id"]),identity)

TTseq.human.refseq.major.isoform.exon.first.non.single.coverage.antisense.corrected = sapply(exon.first.non.single.list,function(x){sum(apply(t(t((t(t(human.refseq.major.isoform.exon.coverage.antisense.corrected.TTseq[x,c("treatment_1.bam","treatment_2.bam")])/alternative.sequencing.depth[c("treatment_1.tabular","treatment_2.tabular")])) 
                                                                                                                                                  - t(alternative.cross.contamination[c("treatment_1.tabular","treatment_2.tabular")]*t((t(t(human.refseq.major.isoform.exon.coverage.antisense.corrected.TTseq[x,c("RNAseq_treatment_1.bam","RNAseq_treatment_2.bam")])/alternative.sequencing.depth[c("RNAseq_treatment_1.tabular","RNAseq_treatment_2.tabular")])))))/(1 - alternative.cross.contamination[c("treatment_1.tabular","treatment_2.tabular")])),1,function(x){x[x < 0] = 0;sum(x)}))})

mNETseq.human.refseq.major.isoform.exon.first.non.single.coverage.antisense.corrected = sapply(exon.first.non.single.list,function(x){sum(human.refseq.major.isoform.exon.coverage.antisense.corrected.mNETseq[x,c("treatment_1.bam","treatment_2.bam")]/size.factors[c("treatment_1.tabular","treatment_2.tabular")])})

elongation.rate.exons.first.non.single = (TTseq.human.refseq.major.isoform.exon.first.non.single.coverage.antisense.corrected/mNETseq.human.refseq.major.isoform.exon.first.non.single.coverage.antisense.corrected)*(mNETseq.size.factor/(10*conversion.factor.to.amount.per.cell))

save(elongation.rate.exons.first.non.single,file="ProcessedData/elongation.rate.exons.first.treatment.RData")

# on first introns (including-single)
human.refseq.major.isoform.intron.first.non.single = human.refseq.major.isoform.intron[which(human.refseq.major.isoform.intron[,"intron_order"] == "first" | human.refseq.major.isoform.intron[,"intron_order"] == "1"),]
human.refseq.major.isoform.intron.first.non.single = human.refseq.major.isoform.intron.first.non.single[which(as.character(human.refseq.major.isoform.intron.first.non.single[,"transcript_id"]) %in% expressed.TR.5.50),]
human.refseq.major.isoform.intron.first.non.single = human.refseq.major.isoform.intron.first.non.single[which(as.character(human.refseq.major.isoform.intron.first.non.single[,"transcript_id"]) %in% protein.coding),]

intron.first.non.single = as.character(human.refseq.major.isoform.intron.first.non.single[,"id"])	
intron.first.non.single.list = tapply(as.character(human.refseq.major.isoform.intron.first.non.single[,"id"]),INDEX = as.character(human.refseq.major.isoform.intron.first.non.single[,"transcript_id"]),identity)

TTseq.human.refseq.major.isoform.intron.first.non.single.coverage.antisense.corrected = sapply(intron.first.non.single.list,function(x){sum(apply(t(t((t(t(human.refseq.major.isoform.intron.coverage.antisense.corrected.TTseq[x,c("treatment_1.bam","treatment_2.bam")])/alternative.sequencing.depth[c("treatment_1.tabular","treatment_2.tabular")])) 
                                                                                                                                                      - t(alternative.cross.contamination[c("treatment_1.tabular","treatment_2.tabular")]*t((t(t(human.refseq.major.isoform.intron.coverage.antisense.corrected.TTseq[x,c("RNAseq_treatment_1.bam","RNAseq_treatment_2.bam")])/alternative.sequencing.depth[c("RNAseq_treatment_1.tabular","RNAseq_treatment_2.tabular")])))))/(1 - alternative.cross.contamination[c("treatment_1.tabular","treatment_2.tabular")])),1,function(x){x[x < 0] = 0;sum(x)}))})

mNETseq.human.refseq.major.isoform.intron.first.non.single.coverage.antisense.corrected = sapply(intron.first.non.single.list,function(x){sum(human.refseq.major.isoform.intron.coverage.antisense.corrected.mNETseq[x,c("treatment_1.bam","treatment_2.bam")]/size.factors[c("treatment_1.tabular","treatment_2.tabular")])})

elongation.rate.introns.first.non.single = (TTseq.human.refseq.major.isoform.intron.first.non.single.coverage.antisense.corrected/mNETseq.human.refseq.major.isoform.intron.first.non.single.coverage.antisense.corrected)*(mNETseq.size.factor/(10*conversion.factor.to.amount.per.cell))

save(elongation.rate.introns.first.non.single,file="ProcessedData/elongation.rate.introns.first.treatment.RData")
