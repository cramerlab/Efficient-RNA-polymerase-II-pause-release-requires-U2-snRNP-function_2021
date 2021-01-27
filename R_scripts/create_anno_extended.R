	load("~/Annotation/human.refseq.anno.RData")
  human.refseq.extended.list = list()
	human.refseq.extended.list[["+"]] = list()
	human.refseq.extended.list[["-"]] = list()
	
	human.refseq.anno = cbind(human.refseq.anno,"exon_order" = NA,"intron_order" = NA,"donor_junction_order" = NA,"acceptor_junction_order" = NA,"UTR_order" = NA,"CDS_order" = NA,"exon_number" = NA,"intron_number" = NA,"donor_junction_number" = NA,"acceptor_junction_number" = NA,"CDS_number" = NA)
	
	registerDoParallel(cores = mc.cores)
	build.human.refseq.extended.list = function(which.chr,which.strand){
		subset.strand.chr.human.refseq.anno = as.matrix(human.refseq.anno[which(human.refseq.anno[,"strand"] == which.strand & human.refseq.anno[,"chr"] == which.chr),])
		human.transcript.ids = unique(as.vector(subset.strand.chr.human.refseq.anno[,"transcript_id"]))
		human.transcript.ids = human.transcript.ids[!is.na(human.transcript.ids)]
		
		final.nr.transcripts = length(human.transcript.ids)
		final.nr.exons = sum(subset.strand.chr.human.refseq.anno[,"type"] == "exon")
		final.nr.CDSs = sum(subset.strand.chr.human.refseq.anno[,"type"] == "CDS")
		
		final.anno = matrix(NA,nrow = final.nr.transcripts + final.nr.exons + 3*(final.nr.exons - final.nr.transcripts) + final.nr.CDSs,ncol = ncol(subset.strand.chr.human.refseq.anno))
		colnames(final.anno) = colnames(subset.strand.chr.human.refseq.anno)
		
		i = 0
		
		for (human.transcript.id in human.transcript.ids){
			subset.transcript.human.refseq.anno = subset.strand.chr.human.refseq.anno[which(subset.strand.chr.human.refseq.anno[,"transcript_id"] %in% human.transcript.id),,drop = FALSE]
			
			subset.exon.human.refseq.anno = subset.transcript.human.refseq.anno[which(subset.transcript.human.refseq.anno[,"type"] == "exon"),,drop = FALSE]
			nr.exons = nrow(subset.exon.human.refseq.anno)
			
			which.cols = c("gene_name","transcript_id")
			
			final.anno[i + 1,c("start","end","type")] = cbind(min(as.numeric(subset.exon.human.refseq.anno[,"start"])),max(as.numeric(subset.exon.human.refseq.anno[,"end"])),"transcript")
			final.anno[i + 1,which.cols] = subset.exon.human.refseq.anno[1,which.cols]
			i = i + 1
			
			if (nr.exons > 1){
				if (which.strand == "+"){
					exon.number = 1:nr.exons
					exon.order = 1:nr.exons
					exon.order[1] = "first"
					exon.order[nr.exons] = "last"
					
					final.anno[i+(1:nr.exons),] = subset.exon.human.refseq.anno
					final.anno[i+(1:nr.exons),"exon_number"] = exon.number[order(as.numeric(subset.exon.human.refseq.anno[,"start"]))]
					final.anno[i+(1:nr.exons),"exon_order"] = exon.order[order(as.numeric(subset.exon.human.refseq.anno[,"start"]))]
					i = i + nr.exons
					
					rest.number = 1:(nr.exons-1)
					rest.order = 1:(nr.exons-1)
					if (nr.exons > 2){
						rest.order[1] = "first"
						rest.order[(nr.exons-1)] = "last"
					}
					
					nr.introns = nr.exons-1
					
					which.cols = c("gene_name","transcript_id")
					
					final.anno[i+(1:nr.introns),c("start","end","type","intron_order","intron_number")] = cbind(as.numeric(subset.exon.human.refseq.anno[,"end"][-nr.exons])+1,as.numeric(subset.exon.human.refseq.anno[,"start"][-1])-1,"intron",rest.order,rest.number)
					final.anno[i+(1:nr.introns)+nr.introns,c("start","end","type","donor_junction_order","donor_junction_number")] = cbind(as.numeric(subset.exon.human.refseq.anno[,"end"][-nr.exons])-1,as.numeric(subset.exon.human.refseq.anno[,"end"][-nr.exons])+2,"donor_junction",rest.order,rest.number)
					final.anno[i+(1:nr.introns)+2*nr.introns,c("start","end","type","acceptor_junction_order","acceptor_junction_number")] = cbind(as.numeric(subset.exon.human.refseq.anno[,"start"][-1])-2,as.numeric(subset.exon.human.refseq.anno[,"start"][-1])+1,"acceptor_junction",rest.order,rest.number)
					final.anno[i+(1:(3*nr.introns)),which.cols] = matrix(rep(subset.exon.human.refseq.anno[1,which.cols],3*nr.introns),nrow = 3*nr.introns,byrow = TRUE)
					i = i + 3*nr.introns
				} else {
					exon.number = 1:nr.exons
					exon.order = 1:nr.exons
					exon.order[1] = "first"
					exon.order[nr.exons] = "last"
					
					subset.exon.human.refseq.anno = subset.exon.human.refseq.anno[order(as.numeric(subset.exon.human.refseq.anno[,"start"])),]
					
					final.anno[i+(1:nr.exons),] = subset.exon.human.refseq.anno
					final.anno[i+(1:nr.exons),"exon_number"] = rev(exon.number[order(as.numeric(subset.exon.human.refseq.anno[,"start"]))])
					final.anno[i+(1:nr.exons),"exon_order"] = rev(exon.order[order(as.numeric(subset.exon.human.refseq.anno[,"start"]))])
					i = i + nr.exons
					
					rest.number = 1:(nr.exons-1)
					rest.order = 1:(nr.exons-1)
					if (nr.exons > 2){
						rest.order[1] = "first"
						rest.order[(nr.exons-1)] = "last"
					}
					
					nr.introns = nr.exons-1
					
					which.cols = c("gene_name","transcript_id")
					
					final.anno[i+(1:nr.introns),c("start","end","type","intron_order","intron_number")] = cbind(as.numeric(subset.exon.human.refseq.anno[,"end"][-nr.exons])+1,as.numeric(subset.exon.human.refseq.anno[,"start"][-1])-1,"intron",rev(rest.order),rev(rest.number))
					final.anno[i+(1:nr.introns)+nr.introns,c("start","end","type","acceptor_junction_order","acceptor_junction_number")] = cbind(as.numeric(subset.exon.human.refseq.anno[,"end"][-nr.exons])-1,as.numeric(subset.exon.human.refseq.anno[,"end"][-nr.exons])+2,"acceptor_junction",rev(rest.order),rev(rest.number))
					final.anno[i+(1:nr.introns)+2*nr.introns,c("start","end","type","donor_junction_order","donor_junction_number")] = cbind(as.numeric(subset.exon.human.refseq.anno[,"start"][-1])-2,as.numeric(subset.exon.human.refseq.anno[,"start"][-1])+1,"donor_junction",rev(rest.order),rev(rest.number))
					final.anno[i+(1:(3*nr.introns)),which.cols] = matrix(rep(subset.exon.human.refseq.anno[1,which.cols],3*nr.introns),nrow = 3*nr.introns,byrow = TRUE)
					i = i + 3*nr.introns
				}
			} else {
				final.anno[i+(1:nr.exons),] = subset.exon.human.refseq.anno
				final.anno[i+(1:nr.exons),"exon_number"] = 1
				final.anno[i+(1:nr.exons),"exon_order"] = "single"
				i = i + nr.exons
			}
			
			subset.CDS.human.refseq.anno = subset.transcript.human.refseq.anno[which(subset.transcript.human.refseq.anno[,"type"] == "CDS"),,drop = FALSE]
			nr.CDSs = nrow(subset.CDS.human.refseq.anno)
			
			if (nr.CDSs > 0){
				if (which.strand == "+"){
					CDS.number = 1:nr.CDSs
					CDS.order = 1:nr.CDSs
					if (nr.CDSs > 1){
						CDS.order[1] = "first"
						CDS.order[nr.CDSs] = "last"
					} else {CDS.order = "single"}
					
					final.anno[i+(1:nr.CDSs),] = subset.CDS.human.refseq.anno
					final.anno[i+(1:nr.CDSs),"CDS_number"] = CDS.number[order(as.numeric(subset.CDS.human.refseq.anno[,"start"]))]
					final.anno[i+(1:nr.CDSs),"CDS_order"] = CDS.order[order(as.numeric(subset.CDS.human.refseq.anno[,"start"]))]
					i = i + nr.CDSs
				} else {
					CDS.number = 1:nr.CDSs
					CDS.order = 1:nr.CDSs
					if (nr.CDSs > 1){
						CDS.order[1] = "first"
						CDS.order[nr.CDSs] = "last"
					} else {CDS.order = "single"}
					
					final.anno[i+(1:nr.CDSs),] = subset.CDS.human.refseq.anno
					final.anno[i+(1:nr.CDSs),"CDS_number"] = rev(CDS.number[order(as.numeric(subset.CDS.human.refseq.anno[,"start"]))])
					final.anno[i+(1:nr.CDSs),"CDS_order"] = rev(CDS.order[order(as.numeric(subset.CDS.human.refseq.anno[,"start"]))])
					i = i + nr.CDSs
				}
			}
		}
		final.anno[,"chr"] = which.chr
		final.anno[,"source"] = "REFSEQ_extended"
		final.anno[,"strand"] = which.strand
		final.anno[,"score"] = "."
		final.anno[,"frame"] = "."
		return(final.anno)
	}
	
	human.refseq.extended.list[["+"]] = foreach(n = human.chrs) %dopar% build.human.refseq.extended.list(n,"+")
	names(human.refseq.extended.list[["+"]]) = human.chrs
	
	human.refseq.extended.list[["-"]] = foreach(n = human.chrs) %dopar% build.human.refseq.extended.list(n,"-")
	names(human.refseq.extended.list[["-"]]) = human.chrs
	
	subset.human.refseq.anno = human.refseq.anno[which(human.refseq.anno[,"type"] != "exon" & human.refseq.anno[,"type"] != "CDS"),]
	human.refseq.extended = rbind(subset.human.refseq.anno,Reduce(rbind,human.refseq.extended.list[["+"]]),Reduce(rbind,human.refseq.extended.list[["-"]]))
	
	rownames(human.refseq.extended) = as.character(1:dim(human.refseq.extended)[1])
	human.refseq.extended = as.data.frame(cbind(human.refseq.extended,"id" = rownames(human.refseq.extended)),stringsAsFactor = FALSE)
	human.refseq.extended[,c("start","end","exon_number","length","placing")] = apply(human.refseq.extended[,c("start","end","exon_number","length","placing")],2,as.numeric)
	human.refseq.extended[,"length"] = human.refseq.extended[,"end"] - human.refseq.extended[,"start"] + 1
	human.refseq.extended[,"placing"] = as.character(1:dim(human.refseq.extended)[1])
	human.refseq.extended$gene_name=sub('\\..*', '', human.refseq.extended$gene_name)
	human.refseq.extended[,"gene_name"] = ensembl.biomart.refseq.2.ensembl.gene.id.mapping[human.refseq.extended[,"gene_name"]]
	dim(human.refseq.extended)
	head(human.refseq.extended)
	
	save(human.refseq.extended,file = "~/Annotation/human.refseq.extended.RData")
