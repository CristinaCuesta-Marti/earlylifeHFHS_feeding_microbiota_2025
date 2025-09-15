library(dada2); packageVersion("dada2")


Sys.setenv(MC_CORES=8) #number of threads



path<-"/data/Run2" #Unzipped fastq file

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz$", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

pdf("plotqual_F_run2_maxee24.pdf")
plotQualityProfile(fnFs[1:8])
dev.off()

pdf("plotqual_R_run2_maxee24.pdf")
plotQualityProfile(fnRs[1:8])
dev.off()

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(17,21), truncLen=c(285,210),    #v1 is 290 245
              maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
class(out) #matrix   
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

pdf("ploterr_F_run2_maxee24.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()


# Apply the core sample inference algorithm (dada) to the dereplicated data


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#Merge Paired Reads#
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
class(seqtab) #matrix
dim(seqtab)
saveRDS(seqtab, "seqtabrun2_24.rds")
