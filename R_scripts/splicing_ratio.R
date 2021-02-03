load("Annotation/human.refseq.major.isoform.exon.RData")
load("Annotation/expressed.TR.5.50.RData")
load("Annotation/protein.coding.RData")

human.refseq.major.isoform.exon=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$transcript_id %in% expressed.TR.5.50),]
human.refseq.major.isoform.exon=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$transcript_id %in% protein.coding),]

###FIRST NON SINGLE Splicing rate####
load("ProcessedData/exon.based.spliced.junction.read.count.list_all.non.single.RData")
treatment_1 <- exon.based.spliced.junction.read.count.list$treatment1.bam
treatment_2 <- exon.based.spliced.junction.read.count.list$treatment2.bam
treatment=cbind(treatment_1[which(treatment_1$id %in% treatment_2$id),], treatment_2[which(treatment_2$id %in% treatment_1$id),12:19])
treatment=cbind(treatment[,1:9], (treatment[,12:19]+treatment[,20:27]))
rownames(treatment)=treatment$id
treatment = treatment[which(treatment$five.prime.total>=30),]

#FIRST EXON >100bp
human.refseq.major.isoform.exon_firstsup100=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$exon_order=="first"),]
human.refseq.major.isoform.exon_firstsup100=human.refseq.major.isoform.exon_firstsup100[which(human.refseq.major.isoform.exon_firstsup100$length>=100),]
firstsup100=human.refseq.major.isoform.exon_firstsup100$transcript_id


treatment_first=treatment[which(treatment$transcript %in% firstsup100),]

control_1 <- exon.based.spliced.junction.read.count.list$control1.bam
control_2 <- exon.based.spliced.junction.read.count.list$control2.bam

control=cbind(control_1[which(control_1$id %in% control_2$id),], control_2[which(control_2$id %in% control_1$id),12:19])
control=cbind(control[,1:9], (control[,12:19]+control[,20:27]))
rownames(control)=control$id

control_first=control[which(control$id %in% treatment_firstsup100$id),]

control_first_5prime.spliced = (control_first["five.prime.spliced"]/control_first["five.prime.total"] + (1-(control_first["five.prime.unspliced"]/control_first["five.prime.total"]))) /2
colnames(control_first_5prime.spliced) = "control_5prime.spliced"
treatment_first_5prime.spliced = (treatment_first["five.prime.spliced"]/treatment_first["five.prime.total"] + (1-(treatment_first["five.prime.unspliced"]/treatment_first["five.prime.total"]))) /2
colnames(treatment_first_5prime.spliced) = "treatment_5prime.spliced"

junctions_first=as.matrix(cbind(control_first_5prime.spliced,treatment_first_5prime.spliced))

###INTERMEDIATE NON SINGLE Splicing rate####
human.refseq.major.isoform.exon_int=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$exon_order!="first" & human.refseq.major.isoform.exon$exon_order!="last"),]
intermediate=human.refseq.major.isoform.exon_int$transcript_id

treatment_int=treatment[which(treatment$transcript %in% intermediate),]
control_int=control[which(control$transcript %in% intermediate),]

control_int_5prime.spliced = (control_int["five.prime.spliced"]/control_int["five.prime.total"] + (1-(control_int["five.prime.unspliced"]/control_int["five.prime.total"]))) /2
colnames(control_int_5prime.spliced) = "control_5prime.spliced"
treatment_int_5prime.spliced = (treatment_int["five.prime.spliced"]/treatment_int["five.prime.total"] + (1-(treatment_int["five.prime.unspliced"]/treatment_int["five.prime.total"]))) /2
colnames(treatment_int_5prime.spliced) = "treatment_5prime.spliced"
control_int_3prime.spliced = (control_int["three.prime.spliced"]/control_int["three.prime.total"] + (1-(control_int["three.prime.unspliced"]/control_int["three.prime.total"]))) /2
colnames(control_int_3prime.spliced) = "control_3prime.spliced"
treatment_int_3prime.spliced = (treatment_int["three.prime.spliced"]/treatment_int["three.prime.total"] + (1-(treatment_int["three.prime.unspliced"]/treatment_int["three.prime.total"]))) /2
colnames(treatment_int_3prime.spliced) = "treatment_3prime.spliced"

junctions_int=as.matrix(cbind(control_int_5prime.spliced,treatment_int_5prime.spliced,control_int_3prime.spliced, treatment_int_3prime.spliced))

###LAST NON SINGLE Splicing rate####
human.refseq.major.isoform.exon_last=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$exon_order=="last"),]
last=human.refseq.major.isoform.exon_last$transcript_id

treatment_last=treatment[which(treatment$transcript %in% last),]
control_last=control[which(control$transcript %in% last),]

control_last_3prime.spliced = (control_last["three.prime.spliced"]/control_last["three.prime.total"] + (1-(control_last["three.prime.unspliced"]/control_last["three.prime.total"]))) /2
colnames(control_last_3prime.spliced) = "control_3prime.spliced"
treatment_last_3prime.spliced = (treatment_last["three.prime.spliced"]/treatment_last["three.prime.total"] + (1-(treatment_last["three.prime.unspliced"]/treatment_last["three.prime.total"]))) /2
colnames(treatment_last_3prime.spliced) = "treatment_3prime.spliced"

junctions_int=as.matrix(cbind(control_last_5prime.spliced,treatment_last_5prime.spliced,control_last_3prime.spliced, treatment_last_3prime.spliced))

#####PLOTING ALL THE JUNCTIONS TOGETHER#########
library(scales)
junctions_all_5prime=rbind(junctions_int[,1:2], junctions_first)
junctions_all_3prime=rbind(junctions_int[,3:4], junctions_last)
dim(junctions_all_5prime)
dim(junctions_all_3prime)

dir.create("Visualization")
pdf(file="Visualization/Junctions_ALL.junctions.pdf")
boxplot(junctions_all_5prime[,"control_5prime.spliced"], junctions_all_5prime[, "treatment_5prime.spliced"], junctions_all_3prime[,"control_3prime.spliced"], junctions_all_3prime[,"treatment_3prime.spliced"], ylim=c(0,1.2),outline=FALSE, notch=TRUE, main= "Splicing rate \n 5'SS 3'SS \n", names=c("control \n 5'SS","treatment \n 5'SS","control \n 3'SS","treatment \n 3'SS"), col=c("darkgrey","indianred3"))
dev.off() 

five.prime=wilcox.test(junctions_all_5prime[,"control_5prime.spliced"], junctions_all_5prime[, "treatment_5prime.spliced"], paired=TRUE, alternative="two.sided")
five.prime
three.prime=wilcox.test(junctions_all_3prime[,"control_3prime.spliced"], junctions_all_3prime[, "treatment_3prime.spliced"], paired=TRUE, alternative="two.sided")
three.prime

library(dplyr)
#####PLOTING THE JUNCTIONS BY EXON ORDER
human.refseq.major.isoform.exon=human.refseq.major.isoform.exon %>%
  group_by(transcript_id) %>%
  mutate(nr_exons = n())

#Transcripts with at least 4 exons
human.refseq.major.isoform.exon_4=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$nr_exons>=4),]
dim(human.refseq.major.isoform.exon_4)
length(unique(human.refseq.major.isoform.exon_4$transcript_id))
tr_4exon=human.refseq.major.isoform.exon_4$id

#FIRST INTRON 5'SS
junctions_first4=junctions_first[which(rownames(junctions_first) %in% tr_4exon),]

#LAST INTRON 3'SS
junctions_last4=junctions_last[which(rownames(junctions_last) %in% tr_4exon),]

#FIRST INTRON 3'SS
exon_2=human.refseq.major.isoform.exon_4[which(human.refseq.major.isoform.exon_4$exon_order==2),]
exon_2=exon_2$id
junctions_2=junctions_int[which(rownames(junctions_int) %in% exon_2),] 
junctions_2_4=junctions_2[which(rownames(junctions_2)%in% tr_4exon),]

#INTERMEDIATE INTRONS
junctions_int4=junctions[which(rownames(junctions_int) %in% tr_4exon),]

#LAST INTRON 5'SS
exon_last_1=human.refseq.major.isoform.exon_4[which(human.refseq.major.isoform.exon_4$exon_order==((human.refseq.major.isoform.exon_4$nr_exons)-1)),]
exon_last_1=exon_last_1$id
length(exon_last_1) 
junctions_last_1=junctions[which(rownames(junctions) %in% exon_last_1),]


pdf(file="Visualization/Junctions_first_int_last.pdf", width=20, height=10)
par(cex.axis=0.8)
boxplot(junctions_first4[,c("control_5prime.spliced")],junctions_first4[,c("treatment_5prime.spliced")],junctions_2_4[,"control_3prime.spliced"],junctions_2_4[,"treatment_3prime.spliced"], junctions_int4[,"control_5prime.spliced"], junctions_int4[,"treatment_5prime.spliced"],junctions_int4[,"control_3prime.spliced"], junctions_int4[,"treatment_3prime.spliced"], junctions_last_1[,"control_5prime.spliced"],junctions_last_1[,"treatment_5prime.spliced"],junctions_last4[,c("control_3prime.spliced")],junctions_last4[,c("treatment_3prime.spliced")], ylim=c(0,1), outline=FALSE, notch=TRUE, main= "Splicing rate \n First, Intermediate, Last \n Ctr vs AMOU2", names=c("control \n First donor","treatment \n First donor","control \n First acceptor", "treatment \n First acceptor","control \n Intermediate donor","treatment \n Intermediate donor", "control \n Intermediate acceptor", "treatment \n Intermediate acceptor","control \n Last donor", "treatment \n Last donor","control \n Last acceptor","treatment \n Last acceptor"), col=c("darkgrey","indianred3", "grey85","rosybrown1"), at = c(1:4, 6:9, 11:14))
dev.off()


five.prime.first=wilcox.test(junctions_first4[,"control_5prime.spliced"], junctions_first4[, "treatment_5prime.spliced"], paired=TRUE, alternative="two.sided")
five.prime.first
three.prime.first=wilcox.test(junctions_2_4[,"control_3prime.spliced"], junctions_2_4[, "treatment_3prime.spliced"], paired=TRUE, alternative="two.sided")
three.prime.first
five.prime.int=wilcox.test(junctions_int4[,"control_5prime.spliced"], junctions_int4[, "treatment_5prime.spliced"], paired=TRUE, alternative="two.sided")
five.prime.int
three.prime.int=wilcox.test(junctions_int4[,"control_3prime.spliced"], junctions_int4[, "treatment_3prime.spliced"], paired=TRUE, alternative="two.sided")
three.prime.int
five.prime.last=wilcox.test(junctions_last_1[,"control_5prime.spliced"], junctions_last_1[, "treatment_5prime.spliced"], paired=TRUE, alternative="two.sided")
five.prime.last
three.prime.last=wilcox.test(junctions_last4[,"control_3prime.spliced"], junctions_last4[, "treatment_3prime.spliced"], paired=TRUE, alternative="two.sided")
three.prime.last
