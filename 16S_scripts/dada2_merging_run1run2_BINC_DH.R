library(dada2); packageVersion("dada2")


Sys.setenv(MC_CORES=8) #number of threads



# Merge multiple runs (if necessary)
st1 <- readRDS("/data/seqtabrun1_24.rds")
st2 <- readRDS("/data/seqtabrun2_24.rds")
st.all <- mergeSequenceTables(st1, st2)


#Remove chimeras

seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# % of Reads Remaining 
sum(seqtab.nochim)/sum(st.all)



##assigning taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "/data/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/data/silva_species_assignment_v138.1.fa.gz")

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

write.csv(seqtab.nochim, file="seqtab_nochim_merged_24.csv")
write.csv(taxa.print, file="taxa_print_merged_24.csv")
write.csv(taxa, file="taxa_merged_24.csv")

