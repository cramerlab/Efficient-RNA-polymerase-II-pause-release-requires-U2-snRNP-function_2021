#STEP1 - Select intronless genes

##SELECT GENES WITH ONE EXON BASED ON NCBI CURATED HG38
load("Annotation/human.refseq.anno.RData")
load("Annotation/human.refseq.major.TR.RData")

#Check how many exons a gene has
intronless_genes_exons=human.refseq.anno[which(human.refseq.anno$type=="exon"),]
intronless_genes=intronless_genes_exons %>%
  group_by(gene_name) %>%
  tally()

#Keep only genes with 1 exon
intronless_genes=intronless_genes[which(intronless_genes$n==1),]
intronless_genes=intronless_genes$gene_name

#Get the intronless genes in our major isoform annotation
intronless_anno=human.refseq.major.TR[which(human.refseq.major.TR$gene_name %in% intronless_genes),]
length(unique(intronless_anno$transcript_id)) 
length(intronless_anno$transcript_id) 
intronless=intronless_anno$transcript_id

#STEP2 - Keep only intronless that do not overlapp with intron-containing annotated genes (1kb)
load("Annotation/human.refseq.extended.RData")
human.refseq.extended.tr=human.refseq.extended[which(human.refseq.extended[,"type"] == "transcript"),]

#All intronless genes (being or not part of the major isoform annotation)
genes_exons=human.refseq.extended[which(human.refseq.extended$type=="exon"),]
all.intronless=genes_exons %>%
  group_by(transcript_id) %>%
  tally()
all.intronless=all.intronless[which(all.intronless$n==1),]
length(all.intronless)
all.intronless=all.intronless$transcript_id
length(all.intronless)

#Genes that are not intronless 
human.refseq.extended.tr.not.intronless=human.refseq.extended.tr[which(human.refseq.extended.tr$transcript_id %notin% all.intronless),]
length(unique(human.refseq.extended.tr.not.intronless$transcript_id))

#Major isoform intronless - to check with overlaps with intron containing genes
major.intronless= human.refseq.extended[which(human.refseq.extended[,"type"] == "transcript" & human.refseq.extended[,"transcript_id"] %in% intronless),]

#Check overlaps between major isoform intronless and intron containing genes
not.intronless = makeGRangesFromDataFrame(human.refseq.extended.tr.not.intronless, keep.extra.columns = TRUE)
intronless = makeGRangesFromDataFrame(major.intronless, keep.extra.columns = TRUE)
#Since I want to make sure they have at least a distance of 1kb, I am going to increase 1kb at each size before overlapping
start(intronless)=start(intronless)-1000
end(intronless)=end(intronless)+1000

hits <- findOverlaps(intronless, not.intronless,type = "any", minoverlap=0, select = "all", ignore.strand=TRUE)
ovpairs <- Pairs(intronless,not.intronless, hits=hits)
pint <- pintersect(ovpairs, ignore.strand=TRUE)

ids_with_overlaps=pint$transcript_id
intronless_no_overlaps=major.intronless[which(major.intronless$transcript_id %notin% ids_with_overlaps),]
dim(intronless_no_overlaps) 

#STEP 3 - Remove > 100bp
for (tr in unique(intronless_no_overlaps$transcript_id)){
  intronless_no_overlaps_utr= intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr),]
  if (intronless_no_overlaps_utr[1,"strand"] =="+"){
    #5'UTR
    # + strand
    intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr),"fiveUTR"]=intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr & intronless_no_overlaps$type =="start_codon"), "start"] - intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr & intronless_no_overlaps$type =="transcript"), "start"]
    #3'UTR
    # + strand
    intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr),"threeUTR"]=intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr & intronless_no_overlaps$type =="transcript"), "end"] - intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr & intronless_no_overlaps$type =="stop_codon"), "end"]
    
  }
  
  if (intronless_no_overlaps_utr[1,"strand"] =="-"){
    #5'UTR
    # - strand
    intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr),"fiveUTR"]=intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr & intronless_no_overlaps$type =="transcript"), "end"] - intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr & intronless_no_overlaps$type =="start_codon"), "end"]
    #3'UTR
    # - strand
    intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr),"threeUTR"]=intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr & intronless_extended$type =="stop_codon"), "start"] - intronless_no_overlaps[which(intronless_no_overlaps$transcript_id == tr & intronless_no_overlaps$type =="transcript"), "start"]
    
  }
}

#Keep only transcripts
intronless_no_overlaps=intronless_no_overlaps[which(intronless_no_overlaps$type == "transcript"),]

#All intronless
summary(intronless_no_overlaps$fiveUTR)
summary(intronless_no_overlaps$threeUTR)

#Remove UTRS >100bp
intronless_anno=intronless_no_overlaps[which(intronless_no_overlaps$fiveUTR <100 & intronless_no_overlaps$threeUTR < 100),]
dim(intronless_anno)

save(intronless_anno, file="Annotation/intronless_anno.RData"))

