#libraries needed
library("dplyr")
library("biomaRt")
`%notin%` <- Negate(`%in%`)

load("Utils/human.chrs.RData")

human.refseq.anno = read.delim(file = "~/Annotation/ncbiRefSeqCurated.gtf",sep="\t",row.names = NULL,header = FALSE,stringsAsFactors = FALSE)
dim(human.refseq.anno)
head(human.refseq.anno)
human.refseq.anno = human.refseq.anno[which(human.refseq.anno[,"V1"] %in% human.chrs),]

colnames(human.refseq.anno)[1:8] = c("chr","source","type","start","end","score","strand","frame")

human.refseq.anno=human.refseq.anno %>% separate(V9, c("gene_id", "transcript_id"), sep=";", extra="drop")

human.refseq.anno$gene_id <- substring(test$gene_id, nchar("gene_id "))
human.refseq.anno$gene_id <- substring(test$gene_id, nchar("transcript_id "))

human.refseq.anno = cbind(human.refseq.anno,"length" = abs(human.refseq.anno[,"start"] - human.refseq.anno[,"end"] + 1),"placing" = as.character(1:dim(human.refseq.anno)[1]))

#ADD ensembl gene name
#remove transcript version because ensembl biomart transcript output doesn't have the version
human.refseq.anno$transcript_id=sub('\\..*', '', human.refseq.anno$transcript_id)


#in this analysis was used the version 96, released on April 2019
ensembl <- useEnsembl(biomart = "ensembl", version=96, dataset = "hsapiens_gene_ensembl")

a=getBM(attributes=c('refseq_mrna', 'refseq_ncrna','external_gene_name'), mart= ensembl)

ensembl.biomart.mapping=a[which(a$refseq_mrna != "" | a$refseq_ncrna != ""),]

ensembl.biomart.mapping[which(ensembl.biomart.mapping[,"refseq_mrna"] != ""),"refseq_ncrna"] = ensembl.biomart.mapping[which(ensembl.biomart.mapping[,"refseq_mrna"] != ""),"refseq_mrna"]

ensembl.biomart.mapping=ensembl.biomart.mapping[,2:3]

colnames(ensembl.biomart.mapping)=c("transcript_id", "gene_name")
ensembl.biomart.mapping=distinct(ensembl.biomart.mapping)

#remove cases where the same transcript_id has two genes
trs = ensembl.biomart.mapping %>% 
group_by(transcript_id) %>% 
summarise(count = n_distinct(gene_name))
trs=trs[which(trs$count >1),]
trs=trs$transcript_id
ensembl.biomart.mapping=ensembl.biomart.mapping[which(ensembl.biomart.mapping$transcript_id %notin% trs),]


mapping=ensembl.biomart.refseq.2.ensembl.gene.id.mapping$gene_name
names(mapping)=ensembl.biomart.refseq.2.ensembl.gene.id.mapping$transcript_id
human.refseq.anno$gene_name=mapping[human.refseq.anno$transcript_id]

save(human.refseq.anno,file = "~/Annotation/human.refseq.anno.RData"))
