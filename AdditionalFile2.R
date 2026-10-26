## Larvae & tank all reads DADA2

setwd("/Volumes/Sarah_SSD/embryo_microbiome/larvae_tank_fastq")  ##set a path##
path <- "/Volumes/Sarah_SSD/embryo_microbiome/larvae_tank_fastq"
list.files(path) #identify files in working directory
library(dada2)

#identify forward and reverse sequence files
fnFs <- sort(list.files(path, pattern="-R1.fastq", full.name = TRUE))
fnRs <- sort(list.files(path, pattern="-R2.fastq", full.names = TRUE))

#extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

#plot quality of sequences
plotQualityProfile(fnFs[1:13]) # visually assess quality 
plotQualityProfile(fnRs[1:13])

#make filtered sub directory
filtFs <-file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <-file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#filter and trim reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),trimLeft=20, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
print(out) #check lost reads

#learn error rates  
errF <-learnErrors(filtFs, multithread = TRUE)
errR <-learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#dereplicate (unique seqs)
derepFs <-derepFastq(filtFs, verbose=TRUE)
derepRs <-derepFastq(filtRs, verbose=TRUE)

names(derepFs) <-sample.names
names(derepRs) <-sample.names

#denoising
dadaFs <-dada(derepFs, err=errF, multithread=TRUE)
dadaRs <-dada(derepRs, err=errR, multithread = TRUE)

dadaFs[[1]]
dadaRs[[1]]

#merge paired reads
mergers <-mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

#make sequence table
seqtab_lar_tank <-makeSequenceTable(mergers)
dim(seqtab_lar_tank)
table(nchar(getSequences(seqtab_lar_tank)))

#remove chimeras
seqtab_LT.nochim <-removeBimeraDenovo(seqtab_lar_tank, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_LT.nochim)
sum(seqtab_LT.nochim)/sum(seqtab_lar_tank)
write.csv(seqtab_LT.nochim, file="larvae_tank_asv.csv") #save ASV file
saveRDS(seqtab_LT.nochim, "seqtab_LT.nochim.rds") #save rds file for merging with embryos

#tracking how sequences went
getN <-function(x) sum(getUniques(x))
track <-cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_LT.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(track)
write.csv(track, file="larvae_tank_tracking.csv")

