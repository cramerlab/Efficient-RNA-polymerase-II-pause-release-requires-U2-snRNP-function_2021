library(GenomicAlignments)

load("Annotation/human.refseq.major.TR.RData")
load("Annotation/human.refseq.major.isoform.exon.RData")
load("Annotation/human.refseq.extended.RData")
load("Utils/human.chrs.RData")

human.refseq.intron.anno.ranges = as(human.refseq.extended[which(human.refseq.extended[,"type"] == "intron"),],"GRanges")

#expressed.TR.5.50 calculation can be found in https://github.com/cramerlab/TT-seq_analysis/
load("ProcessedData/TTseq/expressed.TR.5.50.RData")
load("Annotation/protein.coding.RData")
  
human.refseq.major.isoform.exon=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$transcript_id %in% expressed.TR.5.50),]
human.refseq.major.isoform.exon=human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon$transcript_id %in% protein.coding),]
human.refseq.major.isoform.exon=as.data.frame(human.refseq.major.isoform.exon)
  
human.refseq.major.isoform.exon.all.non.single = human.refseq.major.isoform.exon[which(human.refseq.major.isoform.exon[,"exon_order"] != "single"),]
rownames(human.refseq.major.isoform.exon.all.non.single) = human.refseq.major.isoform.exon.all.non.single[,"id"]
human.refseq.major.isoform.exon.all.non.single = cbind(human.refseq.major.isoform.exon.all.non.single,"id" = human.refseq.major.isoform.exon.all.non.single[,"id"])

exon.based.spliced.junction.read.count.list = list()

bam.files <- dir("TTseq/Bamfiles")

for (bam.file in bam.files){
  build.exon.based.spliced.junction.read.count.list = function(which.chr){
    refseq.constitutive.exons.all.non.single = human.refseq.major.isoform.exon.all.non.single[which(human.refseq.major.isoform.exon.all.non.single[,"chr"] == which.chr),]
    if (dim(refseq.constitutive.exons.all.non.single)[1] > 0){
      junction.counts = cbind(refseq.constitutive.exons.all.non.single,"three.prime.unspliced" = 0,"three.prime.spliced" = 0,"three.prime.total" = 0,"three.prime.inner.mate.exclusive" = 0,"five.prime.unspliced" = 0,"five.prime.spliced" = 0,"five.prime.total" = 0,"five.prime.inner.mate.exclusive" = 0)
      rownames(junction.counts) = rownames(refseq.constitutive.exons.all.non.single)
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = "TTseq/Bamfiles/bam.file",param = param)
      bam = bam[start(left(bam)) <= end(right(bam))]

      bam.single.reads = c(left(bam),right(bam))

      three.prime.junction.ranges = GRanges(seqnames = which.chr,strand = junction.counts[,"strand"],ranges = IRanges(start = junction.counts[,"start"] - 2,end = junction.counts[,"start"] + 1))
      names(three.prime.junction.ranges) = rownames(junction.counts)
      min.2.junction.overlaps = findOverlaps(three.prime.junction.ranges,bam.single.reads,minoverlap = 2L,type = "any",select = "all",ignore.strand = FALSE)
      min.3.junction.overlaps = findOverlaps(three.prime.junction.ranges,bam.single.reads,minoverlap = 3L,type = "any",select = "all",ignore.strand = FALSE)
      min.2.counts = summary(as.factor(names(three.prime.junction.ranges[queryHits(min.2.junction.overlaps)])),maxsum = length(three.prime.junction.ranges))
      min.3.counts = summary(as.factor(names(three.prime.junction.ranges[queryHits(min.3.junction.overlaps)])),maxsum = length(three.prime.junction.ranges))
      junction.counts[names(min.2.counts),"three.prime.total"] = min.2.counts
      junction.counts[names(min.3.counts),"three.prime.unspliced"] = min.3.counts
      junction.counts[,"three.prime.spliced"] = junction.counts[,"three.prime.total"] - junction.counts[,"three.prime.unspliced"]

      five.prime.junction.ranges = GRanges(seqnames = which.chr,strand = junction.counts[,"strand"],ranges = IRanges(start = junction.counts[,"end"] - 1,end = junction.counts[,"end"] + 2))
      names(five.prime.junction.ranges) = rownames(junction.counts)
      min.2.junction.overlaps = findOverlaps(five.prime.junction.ranges,bam.single.reads,minoverlap = 2L,type = "any",select = "all",ignore.strand = FALSE)
      min.3.junction.overlaps = findOverlaps(five.prime.junction.ranges,bam.single.reads,minoverlap = 3L,type = "any",select = "all",ignore.strand = FALSE)
      min.2.counts = summary(as.factor(names(five.prime.junction.ranges[queryHits(min.2.junction.overlaps)])),maxsum = length(five.prime.junction.ranges))
      min.3.counts = summary(as.factor(names(five.prime.junction.ranges[queryHits(min.3.junction.overlaps)])),maxsum = length(five.prime.junction.ranges))
      junction.counts[names(min.2.counts),"five.prime.total"] = min.2.counts
      junction.counts[names(min.3.counts),"five.prime.unspliced"] = min.3.counts
      junction.counts[,"five.prime.spliced"] = junction.counts[,"five.prime.total"] - junction.counts[,"five.prime.unspliced"]

      bam = bam[(start(right(bam)) - end(left(bam))) >= 2]
      inner.mate.granges = GRanges(seqnames = which.chr,strand = strand(bam),ranges = IRanges(start = end(left(bam)) + 1,end = start(right(bam)) - 1))
      refseq.intron.anno.ranges = human.refseq.intron.anno.ranges[seqnames(human.refseq.intron.anno.ranges) == which.chr]
      inner.mate.intron.overlaps = findOverlaps(refseq.intron.anno.ranges,inner.mate.granges,maxgap = 0L,minoverlap = 1L,type = "within",select = "all",ignore.strand = FALSE)
      inner.mate.granges = inner.mate.granges[setdiff(1:length(inner.mate.granges),subjectHits(inner.mate.intron.overlaps))]
      inner.mate.five.prime.overlaps = findOverlaps(five.prime.junction.ranges,inner.mate.granges,minoverlap = 3L,type = "any",select = "all",ignore.strand = FALSE)
      inner.mate.three.prime.overlaps = findOverlaps(three.prime.junction.ranges,inner.mate.granges,minoverlap = 3L,type = "any",select = "all",ignore.strand = FALSE)
      inner.mate.five.prime.counts = summary(as.factor(names(five.prime.junction.ranges[queryHits(inner.mate.five.prime.overlaps)])),maxsum = length(five.prime.junction.ranges))
      inner.mate.three.prime.counts = summary(as.factor(names(three.prime.junction.ranges[queryHits(inner.mate.three.prime.overlaps)])),maxsum = length(three.prime.junction.ranges))
      junction.counts[names(inner.mate.five.prime.counts),"five.prime.inner.mate.exclusive"] = inner.mate.five.prime.counts
      junction.counts[names(inner.mate.three.prime.counts),"three.prime.inner.mate.exclusive"] = inner.mate.three.prime.counts

      junction.counts.backup = junction.counts
      junction.counts[which(junction.counts[,"strand"] == "-"),c("three.prime.unspliced","three.prime.spliced","three.prime.total","three.prime.inner.mate.exclusive")] = junction.counts[which(junction.counts[,"strand"] == "-"),c("five.prime.unspliced","five.prime.spliced","five.prime.total","five.prime.inner.mate.exclusive")]
      junction.counts[which(junction.counts[,"strand"] == "-"),c("five.prime.unspliced","five.prime.spliced","five.prime.total","five.prime.inner.mate.exclusive")] = junction.counts.backup[which(junction.counts[,"strand"] == "-"),c("three.prime.unspliced","three.prime.spliced","three.prime.total","three.prime.inner.mate.exclusive")]

    } else {junction.counts = cbind()}
    return(junction.counts)
  }
  exon.based.spliced.junction.read.count.list.chr = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.exon.based.spliced.junction.read.count.list(n)
  exon.based.spliced.junction.read.count.list[[bam.file]] = do.call("rbind",exon.based.spliced.junction.read.count.list.chr)
}
  
save(exon.based.spliced.junction.read.count.list,file="ProcessedData/TTseq/exon.based.spliced.junction.read.count.list_all.non.single.RData"))
  
