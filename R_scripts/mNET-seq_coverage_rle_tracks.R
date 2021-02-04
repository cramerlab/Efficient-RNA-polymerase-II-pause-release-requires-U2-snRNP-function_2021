library(GenomicAlignments)

load("Utils/human.chrs.RData")
load("Utils/human.chrs.length.RData")


dir.create("ProcessedData/mNETseq")
dir.create("ProcessedData/mNETseq/FragmentEndRleTracks")

bam.files <- dir("mNETseq/Bamfiles")
for (bam.file in bam.files){
		dir.create(file.path("ProcessedData","mNETseq","FragmentEndRleTracks",bam.file))
		build.caizzi.fragment.end.track.list = function(which.chr){
			fragment.end.track.list.chr = list()
			param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
			bam = readGAlignmentPairs(file = file.path("mNETseq","Bamfiles",bam.file),param = param)
			bam = bam[start(left(bam)) <= end(right(bam))]			
			
			end.points = end(right(bam[strand(bam) == "+"]))
			fragment.count = rep(1,length(end.points))
			fragment.counts.aggregate = aggregate(fragment.count,list(id = end.points),function(x){sum(x,na.rm = TRUE)})
			fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
			fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
			fragment.end.track.list.chr[["+"]] = fragment.counts.vec
			
			end.points = start(left(bam[strand(bam) == "-"]))
			fragment.count = rep(1,length(end.points))
			fragment.counts.aggregate = aggregate(fragment.count,list(id = end.points),function(x){sum(x,na.rm = TRUE)})
			fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
			fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
			fragment.end.track.list.chr[["-"]] = fragment.counts.vec
			
			save(fragment.end.track.list.chr,file = file.path("ProcessedData","mNETseq","FragmentEndRleTracks",bam.file,paste0("fragment.end.track.list.",which.chr,".RData")))
			return()
		}
		
		fragment.end.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.fragment.end.track.list(n)
}
	
