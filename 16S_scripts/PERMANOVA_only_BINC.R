library(dada2)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(readxl)
library(tidyverse)
library(ggpubr)
library(vegan)
library(reshape2)
library(tidyr)
library(emmeans)
library(car)
library(agricolae)
library(writexl)
library(broom)
library(svglite)
library(devtools)
library(remotes)
library(permute)
library(pairwiseAdonis)
#install.packages('devtools') # Assuming you don't have devtools package already
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

setwd("C:/Users/DHedayatpour/OneDrive - University College Cork/ADATA onedrive/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items")
samdf <- read_excel("Metadata_BINC_maxee24_ps.xlsx")
sam2 <- samdf %>% remove_rownames %>% column_to_rownames(var="sample")
sam2$sample <- rownames(sam2)

sam2 <- sam2 %>%
  #separate(diet, into = c("diet", "treatment"), sep = " ", extra = "merge") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS", "HFHS B. longum APC1472")),
         #diet = factor(diet, levels = c("Control", "HFHS")),
         #treatment = factor(treatment, levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         sex = factor(sex, levels = c("F", "M")),
         timepoint = factor(timepoint, levels = c("w5", "w10")))

sam2$sex <- factor(sam2$sex, levels = c("M", "F"), labels = c("Male", "Female"))
sam2$timepoint <- factor(sam2$timepoint, levels = c("w5", "w10"), labels = c("week 5", "week 10"))

taxa <- read.csv("taxa_maxee24_ps.csv", sep = ",", row.names = 1)
taxa_mat <- as.matrix(taxa)

seqtab.nochim <- read.csv("seqtab_nochim_maxee24_ps.csv", sep = ",", row.names = 1)
seqtab.nochim_t <- t(seqtab.nochim)
seq_mat <- as.matrix(seqtab.nochim)

ps <- phyloseq(otu_table(seq_mat, taxa_are_rows=FALSE), 
               sample_data(sam2), 
               tax_table(taxa_mat))

ps

##thus removing those for now
ps.noblank <- subset_samples(ps, sample_data(ps)$Sample_or_control !="control")

ps3 <- ps.noblank %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast" &
      Phylum != "Cyanobacteria/Chloroplast"
  )

ps3

##separate the pups and dams, and make do separate ps based on age column

ps_pup <- subset_samples(ps3, sample_data(ps)$age =="pup")
ps_dam <- subset_samples(ps3, sample_data(ps)$age =="dam")

#Make relative abundances
PSr_pup <- transform_sample_counts(ps_pup, function(x) x/sum(x))
PSr_dam <- transform_sample_counts(ps_dam, function(x) x/sum(x))

diet_colours <- c("Control Vehicle" = "#a0a0a4", "HFHS Vehicle" = "#f94040", "HFHS FOS+GOS" = "#addead", "HFHS B. longum APC1472" = "#e1c180")
#Extract components from the original phyloseq object
otu_table <- otu_table(PSr_pup)
tax_table <- tax_table(PSr_pup)
# Extract and modify the sample data to include the reordered diet levels
sample_d <- sample_data(PSr_pup)
sample_d
sample_d <- as.data.frame(sample_d)
sample_d$diet <- factor(sample_d$diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS", "HFHS B. longum APC1472"))
sample_d$timepoint <- factor(sample_d$timepoint, levels = c("week 5", "week 10"))
sample_d$sex <- factor(sample_d$sex, levels = c("Male", "Female"))
PSr_pup_reordered <- phyloseq(otu_table, tax_table, sample_d)
sample_data(PSr_pup_reordered)
class(PSr_pup_reordered)

#Males, Week 5
ps_male_w5 <- subset_samples(PSr_pup_reordered, sex == "Male" & timepoint == "week 5")
# Males, Week 10
ps_male_w10 <- subset_samples(PSr_pup_reordered, sex == "Male" & timepoint == "week 10")
# Females, Week 5
ps_female_w5 <- subset_samples(PSr_pup_reordered, sex == "Female" & timepoint == "week 5")
# Females, Week 10
ps_female_w10 <- subset_samples(PSr_pup_reordered, sex == "Female" & timepoint == "week 10")

##PERMANOVA for all sexes and timepoints
PSr_pup_bray <- phyloseq::distance(PSr_pup_reordered, method = "bray") #measures bray dissimilarity
samples3 <- data.frame(sample_data(PSr_pup_reordered)) #vegan requires sample_data as dataframe
table(samples3$dam_cage)
table(samples3$cohort)
PERMANOVA_pups_all <- adonis2(PSr_pup_bray ~ diet + cohort + dam_cage, data = samples3) #PERMANOVA tests whether the centroids (average distances) of groups (here, the diet groups) are significantly different.

## only if the above adonis result is significant, I delve deeper into the individual diet differences pairwise
diet_levels_pups_all <- unique(samples3$diet)
PERMANOVA_pairwise_pups_all <- data.frame()

for (i in 1:(length(diet_levels_pups_all) - 1)) {
  for (j in (i + 1):length(diet_levels_pups_all)) {
    subset_idx <- samples3$diet %in% c(diet_levels_pups_all[i], diet_levels_pups_all[j])
    dist_sub <- as.dist(as.matrix(PSr_pup_bray)[subset_idx, subset_idx])
    meta_sub <- samples3[subset_idx, ]
    
    res <- adonis2(dist_sub ~ diet + cohort + dam_cage, data = meta_sub, permutations = 999)
    
    PERMANOVA_pairwise_pups_all <- rbind(PERMANOVA_pairwise_pups_all, data.frame(
      Group1 = diet_levels_pups_all[i],
      Group2 = diet_levels_pups_all[j],
      R2 = res$R2[1],
      p = res$`Pr(>F)`[1]
    ))
  }
}
PERMANOVA_pairwise_pups_all$p.adj <- p.adjust(PERMANOVA_pairwise_pups_all$p, method = "fdr")
PERMANOVA_pairwise_pups_all

##PERMANOVA for male w5
PSr_pup_bray_mw5 <- phyloseq::distance(ps_male_w5, method = "bray") #measures bray dissimilarity
samples3_mw5 <- data.frame(sample_data(ps_male_w5)) #vegan requires sample_data as dataframe
table(samples3_mw5$dam_cage)
table(samples3_mw5$cohort)
PERMANOVA_pups_mw5 <- adonis2(PSr_pup_bray_mw5 ~ diet + cohort + dam_cage, data = samples3_mw5) #PERMANOVA tests whether the centroids (average distances) of groups (here, the diet groups) are significantly different.

## only if the above adonis result is significant, I delve deeper into the individual diet differences pairwise
diet_levels_pups_mw5 <- unique(samples3_mw5$diet)
PERMANOVA_pairwise_pups_mw5 <- data.frame()

for (i in 1:(length(diet_levels_pups_mw5) - 1)) {
  for (j in (i + 1):length(diet_levels_pups_mw5)) {
    subset_idx <- samples3_mw5$diet %in% c(diet_levels_pups_mw5[i], diet_levels_pups_mw5[j])
    dist_sub <- as.dist(as.matrix(PSr_pup_bray_mw5)[subset_idx, subset_idx])
    meta_sub <- samples3_mw5[subset_idx, ]
    
    res <- adonis2(dist_sub ~ diet + cohort + dam_cage, data = meta_sub, permutations = 999)
    
    PERMANOVA_pairwise_pups_mw5 <- rbind(PERMANOVA_pairwise_pups_mw5, data.frame(
      Group1 = diet_levels_pups_mw5[i],
      Group2 = diet_levels_pups_mw5[j],
      R2 = res$R2[1],
      p = res$`Pr(>F)`[1]
    ))
  }
}
PERMANOVA_pairwise_pups_mw5$p.adj <- p.adjust(PERMANOVA_pairwise_pups_mw5$p, method = "fdr")
PERMANOVA_pairwise_pups_mw5

##PERMANOVA for female w5
PSr_pup_bray_fw5 <- phyloseq::distance(ps_female_w5, method = "bray") #measures bray dissimilarity
samples3_fw5 <- data.frame(sample_data(ps_female_w5)) #vegan requires sample_data as dataframe
table(samples3_fw5$dam_cage)
table(samples3_fw5$cohort)
PERMANOVA_pups_fw5 <- adonis2(PSr_pup_bray_fw5 ~ diet + cohort + dam_cage, data = samples3_fw5) #PERMANOVA tests whether the centroids (average distances) of groups (here, the diet groups) are significantly different.

## only if the above adonis result is significant, I delve deeper into the individual diet differences pairwise
diet_levels_pups_fw5 <- unique(samples3_fw5$diet)
PERMANOVA_pairwise_pups_fw5 <- data.frame()

for (i in 1:(length(diet_levels_pups_fw5) - 1)) {
  for (j in (i + 1):length(diet_levels_pups_fw5)) {
    subset_idx <- samples3_fw5$diet %in% c(diet_levels_pups_fw5[i], diet_levels_pups_fw5[j])
    dist_sub <- as.dist(as.matrix(PSr_pup_bray_fw5)[subset_idx, subset_idx])
    meta_sub <- samples3_fw5[subset_idx, ]
    
    res <- adonis2(dist_sub ~ diet + cohort + dam_cage, data = meta_sub, permutations = 999)
    
    PERMANOVA_pairwise_pups_fw5 <- rbind(PERMANOVA_pairwise_pups_fw5, data.frame(
      Group1 = diet_levels_pups_fw5[i],
      Group2 = diet_levels_pups_fw5[j],
      R2 = res$R2[1],
      p = res$`Pr(>F)`[1]
    ))
  }
}
PERMANOVA_pairwise_pups_fw5$p.adj <- p.adjust(PERMANOVA_pairwise_pups_fw5$p, method = "fdr")
PERMANOVA_pairwise_pups_fw5

##PERMANOVA for male w10
PSr_pup_bray_mw10 <- phyloseq::distance(ps_male_w10, method = "bray") #measures bray dissimilarity
samples3_mw10 <- data.frame(sample_data(ps_male_w10)) #vegan requires sample_data as dataframe
table(samples3_mw10$dam_cage)
table(samples3_mw10$cohort)
PERMANOVA_pups_mw10 <- adonis2(PSr_pup_bray_mw10 ~ diet + cohort + dam_cage, data = samples3_mw10) #PERMANOVA tests whether the centroids (average distances) of groups (here, the diet groups) are significantly different.

## only if the above adonis result is significant, I delve deeper into the individual diet differences pairwise
diet_levels_pups_mw10 <- unique(samples3_mw10$diet)
PERMANOVA_pairwise_pups_mw10 <- data.frame()

for (i in 1:(length(diet_levels_pups_mw10) - 1)) {
  for (j in (i + 1):length(diet_levels_pups_mw10)) {
    subset_idx <- samples3_mw10$diet %in% c(diet_levels_pups_mw10[i], diet_levels_pups_mw10[j])
    dist_sub <- as.dist(as.matrix(PSr_pup_bray_mw10)[subset_idx, subset_idx])
    meta_sub <- samples3_mw10[subset_idx, ]
    
    res <- adonis2(dist_sub ~ diet + cohort + dam_cage, data = meta_sub, permutations = 999)
    
    PERMANOVA_pairwise_pups_mw10 <- rbind(PERMANOVA_pairwise_pups_mw10, data.frame(
      Group1 = diet_levels_pups_mw10[i],
      Group2 = diet_levels_pups_mw10[j],
      R2 = res$R2[1],
      p = res$`Pr(>F)`[1]
    ))
  }
}
PERMANOVA_pairwise_pups_mw10$p.adj <- p.adjust(PERMANOVA_pairwise_pups_mw10$p, method = "fdr")
PERMANOVA_pairwise_pups_mw10

##PERMANOVA for female w10
PSr_pup_bray_fw10 <- phyloseq::distance(ps_female_w10, method = "bray") #measures bray dissimilarity
samples3_fw10 <- data.frame(sample_data(ps_female_w10)) #vegan requires sample_data as dataframe
table(samples3_fw10$dam_cage)
table(samples3_fw10$cohort)
PERMANOVA_pups_fw10 <- adonis2(PSr_pup_bray_fw10 ~ diet + cohort + dam_cage, data = samples3_fw10) #PERMANOVA tests whether the centroids (average distances) of groups (here, the diet groups) are significantly different.

## only if the above adonis result is significant, I delve deeper into the individual diet differences pairwise
diet_levels_pups_fw10 <- unique(samples3_fw10$diet)
PERMANOVA_pairwise_pups_fw10 <- data.frame()

for (i in 1:(length(diet_levels_pups_fw10) - 1)) {
  for (j in (i + 1):length(diet_levels_pups_fw10)) {
    subset_idx <- samples3_fw10$diet %in% c(diet_levels_pups_fw10[i], diet_levels_pups_fw10[j])
    dist_sub <- as.dist(as.matrix(PSr_pup_bray_fw10)[subset_idx, subset_idx])
    meta_sub <- samples3_fw10[subset_idx, ]
    
    res <- adonis2(dist_sub ~ diet + dam_cage + cohort, data = meta_sub, permutations = 999)
    
    PERMANOVA_pairwise_pups_fw10 <- rbind(PERMANOVA_pairwise_pups_fw10, data.frame(
      Group1 = diet_levels_pups_fw10[i],
      Group2 = diet_levels_pups_fw10[j],
      R2 = res$R2[1],
      p = res$`Pr(>F)`[1]
    ))
  }
}
PERMANOVA_pairwise_pups_fw10$p.adj <- p.adjust(PERMANOVA_pairwise_pups_fw10$p, method = "fdr")
PERMANOVA_pairwise_pups_fw10

### save them all
setwd("C:/Users/DHedayatpour/OneDrive - University College Cork/ADATA onedrive/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/PERMANOVAcagecohort_NEW")

write_xlsx(PERMANOVA_pups_all, "PERMANOVA_pups_all_cagecohort.xlsx")
write_xlsx(PERMANOVA_pairwise_pups_all, "Ppairwise_pups_all_cagecohort.xlsx")

write_xlsx(PERMANOVA_pups_mw5, "PERMANOVA_pups_mw5_cagecohort.xlsx")
write_xlsx(PERMANOVA_pairwise_pups_mw5, "Ppairwise_pups_mw5_cagecohort.xlsx")

write_xlsx(PERMANOVA_pups_fw5, "PERMANOVA_pups_fw5_cagecohort.xlsx")
write_xlsx(PERMANOVA_pairwise_pups_fw5, "Ppairwise_pups_fw5_cagecohort.xlsx")

write_xlsx(PERMANOVA_pups_fw10, "PERMANOVA_pups_fw10_cagecohort.xlsx")
write_xlsx(PERMANOVA_pairwise_pups_fw10, "Ppairwise_pups_fw10_cagecohort.xlsx")

write_xlsx(PERMANOVA_pups_mw10, "PERMANOVA_pups_mw10_cagecohort.xlsx")
write_xlsx(PERMANOVA_pairwise_pups_mw10, "Ppairwise_pups_mw10_cagecohort.xlsx")

################# Dams
diet_colours <- c("Control Vehicle" = "#a0a0a4", "HFHS Vehicle" = "#f94040", "HFHS FOS+GOS" = "#addead", "HFHS B. longum APC1472" = "#e1c180")
#Extract components from the original phyloseq object
otu_table <- otu_table(PSr_dam)
tax_table <- tax_table(PSr_dam)
# Extract and modify the sample data to include the reordered diet levels
sample_d <- sample_data(PSr_dam)
sample_d
sample_d <- as.data.frame(sample_d)
sample_d$diet <- factor(sample_d$diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS", "HFHS B. longum APC1472"))
PSr_dam_reordered <- phyloseq(otu_table, tax_table, sample_d)
sample_data(PSr_dam_reordered)
class(PSr_dam_reordered)


##PERMANOVA for dams
PSr_dam_bray <- phyloseq::distance(PSr_dam_reordered, method = "bray") #measures bray dissimilarity
samples3_dam <- data.frame(sample_data(PSr_dam_reordered)) #vegan requires sample_data as dataframe
#table(samples3_dam$dam_cage)
#table(samples3_dam$cohort)
PERMANOVA_dam <- adonis2(PSr_dam_bray ~ diet, data = samples3_dam) #PERMANOVA tests whether the centroids (average distances) of groups (here, the diet groups) are significantly different.

## only if the above adonis result is significant, I delve deeper into the individual diet differences pairwise
PERMANOVA_pairwise_dams <- pairwise.adonis(phyloseq::distance(PSr_dam_reordered, method = "bray"), sample_data(PSr_dam_reordered)$diet, perm = 999, p.adjust.m = "fdr")

### save both
write_xlsx(PERMANOVA_dam, "PERMANOVA_dam.xlsx")
write_xlsx(PERMANOVA_pairwise_dams, "PERMANOVA_pairwise_dam.xlsx")
