#EXON COVERAGE
load("Annotation/human.refseq.major.isoform.exon.RData"))

#Exon coverage sense
human.refseq.major.isoform.exon.coverage = list()
bam.files <- dir("../TTseq/Bamfiles")
for (bam.file in bam.files){
  index.subsets = split(1:nrow(human.refseq.major.isoform.exon),paste(as.character(human.refseq.major.isoform.exon[,"strand"]),"_",human.refseq.major.isoform.exon[,"chr"],sep = ""))
  coverage.list = list()
  build.coverage.list = function(j){
    from.transcript = strand.chr.human.refseq.major.isoform.exon[j,"start"]
    to.transcript = strand.chr.human.refseq.major.isoform.exon[j,"end"]
    return(sum(as.vector(strand.chr.coverage.from.bam[from.transcript:to.transcript])))
  }

  for (index.subset in names(index.subsets)){
    load(file.path("..","ProcessedData","TTseq","UniqueTranscribedBasesCoverageRleTracks","bam.file",paste0("unique.transcribed.bases.coverage.track.list.",unlist(strsplit(index.subset,split = "_"))[2],".RData")))
    strand.chr.coverage.from.bam = caizzi.cramer.2019.TTseq.K562.PlaB.unique.paired.non.spliced.fragment.coverage.track.list.chr[[unlist(strsplit(index.subset,split = "_"))[1]]]
    strand.chr.human.refseq.major.isoform.exon = human.refseq.major.isoform.exon[index.subsets[[index.subset]],c("start","end","length")]

    registerDoParallel(cores = mc.cores)
    coverage.list[[index.subset]] = foreach(n = 1:nrow(strand.chr.human.refseq.major.isoform.exon),.noexport = setdiff(ls(),c("strand.chr.coverage.from.bam","strand.chr.human.refseq.major.isoform.exon","build.coverage.list"))) %dopar% build.coverage.list(n)
  }
  human.refseq.major.isoform.exon.coverage[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = FALSE)
}	

human.refseq.major.isoform.exon.coverage = sapply(human.refseq.major.isoform.exon.coverage,c)
rownames(human.refseq.major.isoform.exon.coverage) = as.character(human.refseq.major.isoform.exon[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = FALSE)), "id"])
colnames(human.refseq.major.isoform.exon.coverage) = bam.files

save(human.refseq.major.isoform.exon.coverage,file="../ProcessedData/TTseq/human.refseq.major.isoform.exon.coverage.RData")

#Exon coverage antisense
human.refseq.major.isoform.exon.antisense.coverage = list()
for (bam.file in bam.files){
  index.subsets = split(1:nrow(human.refseq.major.isoform.exon),paste(Vectorize(strand.switch)(as.character(human.refseq.major.isoform.exon[,"strand"])),"_",human.refseq.major.isoform.exon[,"chr"],sep = ""))
  coverage.list = list()
  build.coverage.list = function(j){
    from.transcript = strand.chr.human.refseq.major.isoform.exon[j,"start"]
    to.transcript = strand.chr.human.refseq.major.isoform.exon[j,"end"]
    return(sum(as.vector(strand.chr.coverage.from.bam[from.transcript:to.transcript])))
  }

  for (index.subset in names(index.subsets)){
    load(file.path("..","ProcessedData","TTseq","UniqueTranscribedBasesCoverageRleTracks","bam.file",paste0("unique.transcribed.bases.coverage.track.list.",unlist(strsplit(index.subset,split = "_"))[2],".RData")))
    strand.chr.coverage.from.bam = caizzi.cramer.2019.TTseq.K562.PlaB.unique.paired.non.spliced.fragment.coverage.track.list.chr[[unlist(strsplit(index.subset,split = "_"))[1]]]
    strand.chr.human.refseq.major.isoform.exon = human.refseq.major.isoform.exon[index.subsets[[index.subset]],c("start","end","length")]

    registerDoParallel(cores = mc.cores)
    coverage.list[[index.subset]] = foreach(n = 1:nrow(strand.chr.human.refseq.major.isoform.exon),.noexport = setdiff(ls(),c("strand.chr.coverage.from.bam","strand.chr.human.refseq.major.isoform.exon","build.coverage.list"))) %dopar% build.coverage.list(n)
  }
  human.refseq.major.isoform.exon.antisense.coverage[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = FALSE)
}	
human.refseq.major.isoform.exon.antisense.coverage = sapply(human.refseq.major.isoform.exon.antisense.coverage,c)
rownames(human.refseq.major.isoform.exon.antisense.coverage) = as.character(human.refseq.major.isoform.exon[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = FALSE)),"id"])
colnames(human.refseq.major.isoform.exon.antisense.coverage) = bam.files

save(human.refseq.major.isoform.exon.antisense.coverage,file="../ProcessedData/TTseq/human.refseq.major.isoform.exon.antisense.coverage.RData")

#Exon coverage antisense corrected
#antisense.bias.ratio calculation can be found in https://github.com/cramerlab/TT-seq_analysis/
load("../ProcessedData/antisense.bias.ratio.RData"))
  
human.refseq.major.isoform.exon.coverage.antisense.corrected = t(t(human.refseq.major.isoform.exon.coverage - t(t(human.refseq.major.isoform.exon.antisense.coverage)*antisense.bias.ratio[paste(substr(colnames(human.refseq.major.isoform.exon.coverage),1,nchar(colnames(human.refseq.major.isoform.exon.coverage))-4),".tabular",sep = "")]))/(1 - antisense.bias.ratio[paste(substr(colnames(human.refseq.major.isoform.exon.coverage),1,nchar(colnames(human.refseq.major.isoform.exon.coverage))-4),".tabular",sep = "")]^2))
human.refseq.major.isoform.exon.coverage.antisense.corrected[human.refseq.major.isoform.exon.coverage.antisense.corrected < 0] = 0

save(human.refseq.major.isoform.exon.coverage.antisense.corrected,file="../ProcessedData/TTseq/human.refseq.major.isoform.exon.coverage.antisense.corrected.RData")

#INTRON COVERAGE
load("../Annotation/human.refseq.major.isoform.intron.RData"))

#Intron coverage sense
human.refseq.major.isoform.intron.coverage = list()
bam.files <- dir("../TTseq/Bamfiles")
for (bam.file in bam.files){
  index.subsets = split(1:nrow(human.refseq.major.isoform.intron),paste(as.character(human.refseq.major.isoform.intron[,"strand"]),"_",human.refseq.major.isoform.intron[,"chr"],sep = ""))
  coverage.list = list()
  build.coverage.list = function(j){
    from.transcript = strand.chr.human.refseq.major.isoform.intron[j,"start"]
    to.transcript = strand.chr.human.refseq.major.isoform.intron[j,"end"]
    return(sum(as.vector(strand.chr.coverage.from.bam[from.transcript:to.transcript])))
  }
  
  for (index.subset in names(index.subsets)){
    load(file.path("..","ProcessedData","TTseq","UniqueTranscribedBasesCoverageRleTracks","bam.file",paste0("unique.transcribed.bases.coverage.track.list.",unlist(strsplit(index.subset,split = "_"))[2],".RData")))    
    strand.chr.coverage.from.bam = caizzi.cramer.2019.TTseq.K562.PlaB.unique.paired.non.spliced.fragment.coverage.track.list.chr[[unlist(strsplit(index.subset,split = "_"))[1]]]
    strand.chr.human.refseq.major.isoform.intron = human.refseq.major.isoform.intron[index.subsets[[index.subset]],c("start","end","length")]
    
    registerDoParallel(cores = mc.cores)
    coverage.list[[index.subset]] = foreach(n = 1:nrow(strand.chr.human.refseq.major.isoform.intron),.noexport = setdiff(ls(),c("strand.chr.coverage.from.bam","strand.chr.human.refseq.major.isoform.intron","build.coverage.list"))) %dopar% build.coverage.list(n)
  }
  human.refseq.major.isoform.intron.coverage[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = FALSE)
}	

human.refseq.major.isoform.intron.coverage = sapply(human.refseq.major.isoform.intron.coverage,c)
rownames(human.refseq.major.isoform.intron.coverage) = as.character(human.refseq.major.isoform.intron[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = FALSE)), "id"])
colnames(human.refseq.major.isoform.intron.coverage) = bam.files

save(human.refseq.major.isoform.intron.coverage,file="../ProcessedData/TTseq/human.refseq.major.isoform.intron.coverage.RData")

#Intron coverage antisense
human.refseq.major.isoform.intron.antisense.coverage = list()
for (bam.file in bam.files){
  index.subsets = split(1:nrow(human.refseq.major.isoform.intron),paste(Vectorize(strand.switch)(as.character(human.refseq.major.isoform.intron[,"strand"])),"_",human.refseq.major.isoform.intron[,"chr"],sep = ""))
  coverage.list = list()
  build.coverage.list = function(j){
    from.transcript = strand.chr.human.refseq.major.isoform.intron[j,"start"]
    to.transcript = strand.chr.human.refseq.major.isoform.intron[j,"end"]
    return(sum(as.vector(strand.chr.coverage.from.bam[from.transcript:to.transcript])))
  }
  
  for (index.subset in names(index.subsets)){
    load(load(file.path("../ProcessedData","TTseq","UniqueTranscribedBasesCoverageRleTracks","bam.file",paste0("unique.transcribed.bases.coverage.track.list.",unlist(strsplit(index.subset,split = "_"))[2],".RData")))
    
    strand.chr.coverage.from.bam = caizzi.cramer.2019.TTseq.K562.PlaB.unique.paired.non.spliced.fragment.coverage.track.list.chr[[unlist(strsplit(index.subset,split = "_"))[1]]]
    strand.chr.human.refseq.major.isoform.intron = human.refseq.major.isoform.intron[index.subsets[[index.subset]],c("start","end","length")]
    
    registerDoParallel(cores = mc.cores)
    coverage.list[[index.subset]] = foreach(n = 1:nrow(strand.chr.human.refseq.major.isoform.intron),.noexport = setdiff(ls(),c("strand.chr.coverage.from.bam","strand.chr.human.refseq.major.isoform.intron","build.coverage.list"))) %dopar% build.coverage.list(n)
  }
  human.refseq.major.isoform.intron.antisense.coverage[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = FALSE)
}	
human.refseq.major.isoform.intron.antisense.coverage = sapply(human.refseq.major.isoform.intron.antisense.coverage,c)
rownames(human.refseq.major.isoform.intron.antisense.coverage) = as.character(human.refseq.major.isoform.intron[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = FALSE)),"id"])
colnames(human.refseq.major.isoform.intron.antisense.coverage) = bam.files

save(human.refseq.major.isoform.intron.antisense.coverage,file="../ProcessedData/TTseq/human.refseq.major.isoform.intron.antisense.coverage.RData")

#Intron coverage antisense corrected
#antisense.bias.ratio calculation can be found in https://github.com/cramerlab/TT-seq_analysis/
load("ProcessedData/TTseq/antisense.bias.ratio.RData"))

human.refseq.major.isoform.intron.coverage.antisense.corrected = t(t(human.refseq.major.isoform.intron.coverage - t(t(human.refseq.major.isoform.intron.antisense.coverage)*antisense.bias.ratio[paste(substr(colnames(human.refseq.major.isoform.intron.coverage),1,nchar(colnames(human.refseq.major.isoform.intron.coverage))-4),".tabular",sep = "")]))/(1 - antisense.bias.ratio[paste(substr(colnames(human.refseq.major.isoform.intron.coverage),1,nchar(colnames(human.refseq.major.isoform.intron.coverage))-4),".tabular",sep = "")]^2))
human.refseq.major.isoform.intron.coverage.antisense.corrected[human.refseq.major.isoform.intron.coverage.antisense.corrected < 0] = 0

save(human.refseq.major.isoform.intron.coverage.antisense.corrected,file="../ProcessedData/TTseq/human.refseq.major.isoform.intron.coverage.antisense.corrected.RData")

#SPIKEIN COVERAGE
load(file.path("../AnnotationObjects","spikein.anno.RData"))

spikein.anno.transcript.coverage = list()
for (bam.file in bam.files){
  index.subsets = split(1:nrow(spikein.anno),paste(as.character(spikein.anno[,"strand"]),"_",spikein.anno[,"chr"],sep = ""))
  
  coverage.list = list()
  
  build.coverage.list = function(j){
    from.transcript = strand.chr.spikein.anno[j,"start"]
    to.transcript = strand.chr.spikein.anno[j,"end"]
    return(sum(as.vector(strand.chr.coverage.from.bam[from.transcript:to.transcript])))
  }
  
  for (index.subset in names(index.subsets)){
    load(file.path("..","ProcessedData","TTseq","UniquePairedFragmentSpikeinCoverageRleTracks",bam.file,paste0("unique.paired.fragment.spikein.coverage.track.list.",unlist(strsplit(index.subset,split = "_"))[2],".RData")))
    
    strand.chr.coverage.from.bam = unique.paired.fragment.spikein.coverage.track.list.chr[[unlist(strsplit(index.subset,split = "_"))[1]]]
    strand.chr.spikein.anno = spikein.anno[index.subsets[[index.subset]],c("start","end","length")]
    
    registerDoParallel(cores = mc.cores)
    coverage.list[[index.subset]] = foreach(n = 1:nrow(strand.chr.spikein.anno),.noexport = setdiff(ls(),c("strand.chr.coverage.from.bam","strand.chr.spikein.anno","build.coverage.list"))) %dopar% build.coverage.list(n)
  }
  spikein.anno.transcript.coverage[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = FALSE)
}
spikein.anno.transcript.coverage = sapply(spikein.anno.transcript.coverage,c)
rownames(spikein.anno.transcript.coverage) = as.character(spikein.anno[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = FALSE)),"chr"])
colnames(spikein.anno.transcript.coverage) = bam.files

save(spikein.anno.transcript.coverage,file="../ProcessedData/TTseq/spikein.anno.transcript.coverage.RData")




