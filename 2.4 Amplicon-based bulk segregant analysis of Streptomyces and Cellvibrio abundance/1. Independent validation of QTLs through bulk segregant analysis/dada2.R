source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
library("dada2")

source("http://www.bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")

source("https://bioconductor.org/biocLite.R")
biocLite()

source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 
library(bioconductor); packageVersion("bioconductor")

install.packages("BiocInstaller", repos="http://bioconductor.org/packages/3.3/bioc")
library(BiocInstaller)
biocLite(c("IRanges"))
library(devtools)
devtools::install_github("benjjneb/dada2")

###bacteria####

path <- "/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/B2R-project/victorc/RIL_final_analysis/16s"
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[3:4])



filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,truncLen=c(300,280), trimLeft =c(17,21),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


head(out)
print(out)
#write.csv(out, "out300220.csv")
keep <- out[,"reads.out"] > 100 # Or other cutoff
filtFs <- file.path(fnFs, filtFs)[keep]
filtRs <- file.path(fnRs, filtRs)[keep]

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]




errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:
plotErrors(errF, nominalQ=TRUE)


derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#We are now ready to apply the core sample inference algorithm to the dereplicated data.

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]


#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


#Track reads through the pipeline


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "track2.2.csv")

#Assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)


#species level
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")


write.csv(taxa, "taxa.csv")
write.csv(t(seqtab.nochim), "counts.csv")
write.csv(cbind(t(seqtab.nochim),taxa), "seqtab.csv", quote=FALSE)





###Fungi####
path <- "/mnt/nfs/bioinfdata/home/NIOO/victorc/RIL/RIL_merge/fungi/dada2/"
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])



filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,truncLen=c(300,300), trimLeft =c(18,18),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
print(out)
#keep <- out[,"reads.out"] > 100 # Or other cutoff
#filtFs <- file.path(fnFs, filtFs)[keep]
#filtRs <- file.path(fnRs, filtRs)[keep]




errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:
plotErrors(errF, nominalQ=TRUE)


derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#We are now ready to apply the core sample inference algorithm to the dereplicated data.

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]


#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


#Track reads through the pipeline


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "track2.2.csv")

#Assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "/mnt/nfs/bioinfdata/home/NIOO/victorc/RIL/RIL_merge/fungi/dada2/silva_nr_v132_train_set.fa.gz", multithread=TRUE)



taxa <- addSpecies(taxa, "/mnt/nfs/bioinfdata/home/NIOO/victorc/RIL/RIL_merge/fungi/dada2/silva_species_assignment_v132.fa.gz")


write.csv(taxa, "taxa.csv")

write.csv(t(seqtab.nochim), "counts.csv")

write.csv(cbind(t(seqtab.nochim),taxa), "seqtab.csv", quote=FALSE)

