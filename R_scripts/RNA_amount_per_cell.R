#Calculate spikein.number.per.cell with MW [g/mol]

# Nucleic Acid Molecular Weight Conversions [g/mol]

# Exact M.W. of ssRNA (e.g., RNA Transcript):
# M.W. = (An x 329.2) + (Un x 306.2) + (Cn x 305.2) + (Gn x 345.2) + 159
# An, Un, Cn, and Gn are the number of each respective nucleotide within the polynucleotide.
# M.W. calculated is valid at physiological pH.
# Addition of "159" to the M.W. takes into account the M.W. of a 5' triphosphate.

# Molecular Weight of 4-Thiouridine (Ribonucleotide Monophosphate)
# U + S - O: 306.2 + 32.06 - 16 = 322.26

load("Annotation/spikein.anno.corrected.RData")
spikein.genome = readDNAStringSet(filepath = file.path(prewd,"RawData","SpikeIns","spikeins.corrected.fa"))

nucleotide.counts = letterFrequency(spikein.genome,multi.mer.list[["1"]][1:4])
rownames(nucleotide.counts) = names(spikein.genome)

spikein.molecular.weight = c(nucleotide.counts[c("chrS2","chrS4","chrS8"),"A"]*329.2 + 0.9*nucleotide.counts[c("chrS2","chrS4","chrS8"),"T"]*306.2 + nucleotide.counts[c("chrS2","chrS4","chrS8"),"C"]*305.2 + nucleotide.counts[c("chrS2","chrS4","chrS8"),"G"]*345.2 + 0.1*nucleotide.counts[c("chrS2","chrS4","chrS8"),"T"]*322.26 + 159,nucleotide.counts[c("chrS5","chrS9","chrS12"),"A"]*329.2 + nucleotide.counts[c("chrS5","chrS9","chrS12"),"T"]*306.2 + nucleotide.counts[c("chrS5","chrS9","chrS12"),"C"]*305.2 + nucleotide.counts[c("chrS5","chrS9","chrS12"),"G"]*345.2 + 159)

# labeled spikeins: Spike2, Spike4  & Spike8 (10%)
# 50 x 6 ng spike ins in 5 x 10^7 cells (for all spikeins)
# 50 x 3 ng spike ins in 5 x 10^7 cells (for labeled spikeins)

labeled.spikein.molecular.weight = sum(spikein.molecular.weight[c("chrS2","chrS4","chrS8")])

# amount of substance  = (1ng = 10^(-9) [g]) * (Avogadro constant = 6.02214085774*10^23 [mol^(-1)]) / (molecular weight [g/mol])

spikein.number = 50*10^(-9)*(6.02214085774*10^23)/spikein.molecular.weight[c("chrS2","chrS4","chrS8")] #LIVIA

spikein.number.per.cell = spikein.number/(5.0*10^7)

load("ProcessedData/spikein.anno.transcript.coverage.RData")

conversion.factor.to.amount.per.cell = mean(apply(spikein.anno.transcript.coverage[c("chrS2","chrS4","chrS8"),c("control_1.bam","control_2.bam")]/(spikein.anno.corrected[which(spikein.anno.corrected[,"chr"] %in% c("chrS2","chrS4","chrS8")),"end"]*spikein.number.per.cell),2,median))
conversion.factor.to.amount.per.cell

cv = sd(apply(spikein.anno.transcript.coverage[c("chrS2","chrS4","chrS8"),]/(spikein.anno.corrected[which(spikein.anno.corrected[,"chr"] %in% c("chrS2","chrS4","chrS8")),"end"]*spikein.number.per.cell),2,median))/conversion.factor.to.amount.per.cell

save(conversion.factor.to.amount.per.cell,file = "ProcessedData/conversion.factor.to.amount.per.cell.RData")


