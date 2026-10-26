### R Code for Climate Change Lobster Embryo Project ###
### DADA2 for 32 Embryo Samples, plus a positive and negative sequencing control ###

library(dada2)

rm(list=ls()) ## clean up previous run ##
setwd("/Volumes/Sarah_SSD/embryo_microbiome/embryo_fastq/fastq")  ##set a path
path <- "/Volumes/Sarah_SSD/embryo_microbiome/embryo_fastq/fastq"
list.files(path) #identify files in working directory

#identify forward and reverse sequence files
fnFs <- sort(list.files(path, pattern="-R1.fastq", full.name = TRUE))
fnRs <- sort(list.files(path, pattern="-R2.fastq", full.names = TRUE))

#extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

#plot quality of sequences
plotQualityProfile(fnFs[1:34]) # 1:number of samples
plotQualityProfile(fnRs[1:34])

#make filtered subdirectory
filtFs <-file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <-file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#filter and trim reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,160),trimLeft=20, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
print(out) #check output to ensure not too many sequences lost. 

#save as embryo read table
out_emb <- as.table(out)


#learn error rates  
errF <-learnErrors(filtFs, multithread = TRUE)
errR <-learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#dereplicate
derepFs <-derepFastq(filtFs, verbose=TRUE)
derepRs <-derepFastq(filtRs, verbose=TRUE)

#rename 
names(derepFs) <-sample.names
names(derepRs) <-sample.names

#denoise
dadaFs <-dada(derepFs, err=errF, multithread=TRUE)
dadaRs <-dada(derepRs, err=errR, multithread = TRUE)

dadaFs[[1]]
dadaRs[[1]]

#merge paired reads
mergers <-mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

#make sequence table
seqtab_emb <-makeSequenceTable(mergers)
dim(seqtab_emb)
table(nchar(getSequences(seqtab_emb)))

#remove chimeras
seqtab_emb.nochim <-removeBimeraDenovo(seqtab_emb, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_emb.nochim)
sum(seqtab_emb.nochim)/sum(seqtab_emb)
write.csv(seqtab_emb.nochim, file="emb_asv.csv") #save ASV file
saveRDS(seqtab_emb.nochim, "/Volumes/Sarah_SSD/embryo_microbiome/embryo_fastq/seqtab_emb.nochim.rds") #save rds file for merging with larvae and tank

#tracking how sequences went
getN <-function(x) sum(getUniques(x))
track <-cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_emb.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file="/Volumes/Sarah_SSD/embryo_microbiome/embryo_fastq/embryo_tracking.csv")

