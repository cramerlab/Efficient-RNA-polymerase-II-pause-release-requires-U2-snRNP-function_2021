library(GenomicAlignments)

load("Utils/human.chrs.RData")
load("Utils/human.chrs.length.RData")

dir.create("ProcessedData/TTseq")
dir.create("ProcessedData/TTseq/UniqueTranscribedBasesCoverageRleTracks")

bam.files <- dir("TTseq/Bamfiles")

for (bam.file in bam.files){
  dir.create(file.path("ProcessedData","TTseq","UniqueTranscribedBasesCoverageRleTracks",bam.file)
  build.unique.transcribed.bases.coverage.track.list = function(which.chr){
    unique.transcribed.bases.coverage.track.list.chr = list()
    param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
    bam = readGAlignmentPairs(file = file.path("TTseq","Bamfiles",bam.file),param = param)
    bam = bam[start(left(bam)) <= end(right(bam))]
    rle.vec = Rle(0,human.chrs.lengths[which.chr])
    coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(bam[strand(bam) == "+"])),end = end(right(bam[strand(bam) == "+"])))))[[which.chr]]
    rle.vec[1:length(coverage.vec)] = coverage.vec
    unique.transcribed.bases.coverage.track.list.chr[["+"]] = rle.vec
    rle.vec = Rle(0,human.chrs.lengths[which.chr])
    coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(bam[strand(bam) == "-"])),end = end(right(bam[strand(bam) == "-"])))))[[which.chr]]
    rle.vec[1:length(coverage.vec)] = coverage.vec
    unique.transcribed.bases.coverage.track.list.chr[["-"]] = rle.vec
    save(unique.transcribed.bases.coverage.track.list.chr,file = file.path("ProcessedData","TTseq","UniqueTranscribedBasesCoverageRleTracks","bam.file",paste0("unique.transcribed.bases.coverage.track.list.",which.chr,".RData"))
    return()
  }
    unique.transcribed.bases.coverage.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.unique.transcribed.bases.coverage.track.list(n)
}
  
