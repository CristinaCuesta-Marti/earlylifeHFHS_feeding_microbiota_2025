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

###Making the phyloseq object for downstream analysis###
setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items")

samdf <- read_excel("Metadata_BINC_maxee24_ps.xlsx")

sam2 <- samdf %>% remove_rownames %>% column_to_rownames(var="sample")
sam2$sample <- rownames(sam2)

##load the tax table
taxa <- read.csv("taxa_maxee24_ps.csv", sep = ",", row.names = 1)
taxa_mat <- as.matrix(taxa)

##load the otu table
seqtab.nochim <- read.csv("seqtab_nochim_maxee24_ps.csv", sep = ",", row.names = 1)
seqtab.nochim_t <- t(seqtab.nochim)

seq_mat <- as.matrix(seqtab.nochim)

###make a phyloseq object 

##to make a phyloseq object - we need metadata as df, and other 2 as matrix

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

####################### pups

#############First calculate shannon diversity and assign it to a column in metadata
div.all.1_pup <- estimate_richness(ps_pup, measures = "Shannon")
metadata_pup$shan <- div.all.1_pup

rownames(metadata_pup) == rownames(div.all.1_pup)

##Summary lm-FOS+GOS
mat_restored_FOSGOS  = cbind(c(1/2, -1, 1/2))
rownames(mat_restored_FOSGOS) <- c("Control Vehicle","HFHS Vehicle","HFHS FOS+GOS")
colnames(mat_restored_FOSGOS) <- c("Restored")

lm_summary_pups_FOSGOS <- metadata_pup %>%
  filter(diet != "HFHS B. longum APC1472") %>%
  #mutate(Sample = factor(Sample)) %>%
  #mutate(Phylum = factor(Phylum)) %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS"))) %>%
  #group_by() %>%
  reframe(
    lm(shan$Shannon ~  diet * sex * timepoint, 
       contrasts = list(diet = mat_restored_FOSGOS),
       data = . 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  #ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_summary_pups_FOSGOS
setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Stats")
write_xlsx(lm_summary_pups_FOSGOS, "lm_summary_pups_FOSGOS.xlsx")

##Summary lm-FOS+GOS w5
mat_restored_FOSGOS  = cbind(c(1/2, -1, 1/2))
rownames(mat_restored_FOSGOS) <- c("Control Vehicle","HFHS Vehicle","HFHS FOS+GOS")
colnames(mat_restored_FOSGOS) <- c("Restored")

lm_summary_pups_FOSGOS_w5 <- metadata_pup %>%
  filter(diet != "HFHS B. longum APC1472") %>%
  filter(timepoint != "w10") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS"))) %>%
  #group_by()
  reframe(
    lm(shan$Shannon ~  diet * sex, 
       contrasts = list(diet = mat_restored_FOSGOS),
       data = . 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  #ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_summary_pups_FOSGOS_w5
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Stats")
write_xlsx(lm_summary_pups_FOSGOS_w5, "lm_summary_pups_FOSGOS_w5.xlsx")

##Summary lm-FOS+GOS w10
lm_summary_pups_FOSGOS_w10 <- metadata_pup %>%
  filter(diet != "HFHS B. longum APC1472") %>%
  filter(timepoint != "w5") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS"))) %>%
  #group_by()
  reframe(
    lm(shan$Shannon ~  diet * sex, 
       contrasts = list(diet = mat_restored_FOSGOS),
       data = . 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  #ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_summary_pups_FOSGOS_w10
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Stats")
write_xlsx(lm_summary_pups_FOSGOS_w10, "lm_summary_pups_FOSGOS_w10.xlsx")

#Contrasts lm-FOS+GOS w5 grouped by sex
mat_restored_FOSGOS  = cbind(c(1/2, -1, 1/2))
rownames(mat_restored_FOSGOS) <-  c("Control Vehicle","HFHS Vehicle","HFHS FOS+GOS")
colnames(mat_restored_FOSGOS) <- c("Restored")

lm_contrasts_pups_FOSGOS_w5 <- metadata_pup %>%
  filter(diet != "HFHS B. longum APC1472") %>%
  filter(timepoint != "w10") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS"))) %>%
  group_by(sex) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_FOSGOS),
       data = across(everything()) 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_contrasts_pups_FOSGOS_w5
write_xlsx(lm_contrasts_pups_FOSGOS_w5, "lm_contrasts_pups_FOSGOS_w5.xlsx")

#Contrasts lm-FOS+GOS w10 grouped by sex
mat_restored_FOSGOS  = cbind(c(1/2, -1, 1/2))
rownames(mat_restored_FOSGOS) <-  c("Control Vehicle","HFHS Vehicle","HFHS FOS+GOS")
colnames(mat_restored_FOSGOS) <- c("Restored")

lm_contrasts_pups_FOSGOS_w10 <- metadata_pup %>%
  filter(diet != "HFHS B. longum APC1472") %>%
  filter(timepoint != "w5") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS"))) %>%
  group_by(sex) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_FOSGOS),
       data = across(everything()) 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_contrasts_pups_FOSGOS_w10
write_xlsx(lm_contrasts_pups_FOSGOS_w10, "lm_contrasts_pups_FOSGOS_w10.xlsx")

###Contrasts lm-FOS+GOS grouped by sex
lm_contrasts_pups_FOSGOS_sex <- metadata_pup %>%
  filter(diet != "HFHS B. longum APC1472") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS"))) %>%
  group_by(sex) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_FOSGOS),
       data = across(everything()) 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_contrasts_pups_FOSGOS_sex
write_xlsx(lm_contrasts_pups_FOSGOS_sex, "lm_contrasts_pups_FOSGOS_sex.xlsx")

#Contrasts lm-FOS+GOS grouped by timepoint
mat_restored_FOSGOS  = cbind(c(1/2, -1, 1/2))
rownames(mat_restored_FOSGOS) <-  c("Control Vehicle","HFHS Vehicle","HFHS FOS+GOS")
colnames(mat_restored_FOSGOS) <- c("Restored")

lm_contrasts_pups_FOSGOS_timepoint <- metadata_pup %>%
  filter(diet != "HFHS B. longum APC1472") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS"))) %>%
  group_by(timepoint) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_FOSGOS),
       data = across(everything()) 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_contrasts_pups_FOSGOS_timepoint
write_xlsx(lm_contrasts_pups_FOSGOS_timepoint, "lm_contrasts_pups_FOSGOS_timepoint.xlsx")

##Summary lm 1472
mat_restored_1472  = cbind(c(1/2, -1, 1/2))
rownames(mat_restored_1472) <- c("Control Vehicle","HFHS Vehicle","HFHS B. longum APC1472")
colnames(mat_restored_1472) <- c("Restored")

lm_summary_pups_1472 <- metadata_pup %>%
  filter(diet != "HFHS FOS+GOS") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS B. longum APC1472"))) %>%
  #group_by() %>%
  reframe(
    lm(shan$Shannon ~  diet * sex * timepoint, 
       contrasts = list(diet = mat_restored_1472),
       data = . 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  #ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_summary_pups_1472
write_xlsx(lm_summary_pups_1472, "lm_summary_pups_1472.xlsx")

##Summary lm 1472 w5
mat_restored_1472  = cbind(c(1/2, -1, 1/2))
rownames(mat_restored_1472) <- c("Control Vehicle","HFHS Vehicle","HFHS B. longum APC1472")
colnames(mat_restored_1472) <- c("Restored")

lm_summary_pups_1472_w5 <- metadata_pup %>%
  filter(diet != "HFHS FOS+GOS") %>%
  filter(timepoint != "w10") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS B. longum APC1472"))) %>%
  #group_by() %>%
  reframe(
    lm(shan$Shannon ~  diet * sex, 
       contrasts = list(diet = mat_restored_1472),
       data = . 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  #ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_summary_pups_1472_w5
write_xlsx(lm_summary_pups_1472_w5, "lm_summary_pups_1472_w5.xlsx")

##Summary lm 1472 w10
lm_summary_pups_1472_w10 <- metadata_pup %>%
  filter(diet != "HFHS FOS+GOS") %>%
  filter(timepoint != "w5") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS B. longum APC1472"))) %>%
  #group_by() %>%
  reframe(
    lm(shan$Shannon ~  diet * sex, 
       contrasts = list(diet = mat_restored_1472),
       data = . 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  #ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_summary_pups_1472_w10
write_xlsx(lm_summary_pups_1472_w10, "lm_summary_pups_1472_w10.xlsx")

##Contrasts lm-1472 w5 grouped by sex
lm_contrasts_pups_1472_w5 <- metadata_pup %>%
  filter(diet != "HFHS FOS+GOS") %>%
  filter(timepoint != "w10") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS B. longum APC1472"))) %>%
  group_by(sex) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_1472),
       data = across(everything()) 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_contrasts_pups_1472_w5
write_xlsx(lm_contrasts_pups_1472_w5, "lm_contrasts_pups_1472_w5.xlsx")

##Contrasts lm-1472 w10 grouped by sex
lm_contrasts_pups_1472_w10 <- metadata_pup %>%
  filter(diet != "HFHS FOS+GOS") %>%
  filter(timepoint != "w5") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS B. longum APC1472"))) %>%
  group_by(sex) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_1472),
       data = across(everything()) 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_contrasts_pups_1472_w10
write_xlsx(lm_contrasts_pups_1472_w10, "lm_contrasts_pups_1472_w10.xlsx")

##Contrasts lm-1472 grouped by sex
lm_contrasts_pups_1472_sex <- metadata_pup %>%
  filter(diet != "HFHS FOS+GOS") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS B. longum APC1472"))) %>%
  group_by(sex) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_1472),
       data = across(everything()) 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_contrasts_pups_1472_sex
write_xlsx(lm_contrasts_pups_1472_sex, "lm_contrasts_pups_1472_sex.xlsx")

##Contrasts lm-1472 grouped by timepoint
lm_contrasts_pups_1472_timepoint <- metadata_pup %>%
  filter(diet != "HFHS FOS+GOS") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS B. longum APC1472"))) %>%
  group_by(timepoint) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_1472),
       data = across(everything()) 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_contrasts_pups_1472_timepoint
write_xlsx(lm_contrasts_pups_1472_timepoint, "lm_contrasts_pups_1472_timepoint.xlsx")

##Tukeys pups FOS+GOS grouped by sex
tukey_pups_FOSGOS_sex = metadata_pup %>% 
  filter(diet != "HFHS B. longum APC1472") %>% 
  #mutate(Mouse_ID = factor(Mouse_ID)) %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS FOS+GOS"))) %>% 
  group_by(sex) %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_pups_FOSGOS_sex
write_xlsx(tukey_pups_FOSGOS_sex, "tukey_pups_FOSGOS_sex.xlsx")

##Tukeys pups FOS+GOS grouped by sex w5 only
tukey_pups_FOSGOS_sex_w5 = metadata_pup %>% 
  filter(diet != "HFHS B. longum APC1472") %>% 
  filter(timepoint != "w10") %>%
  #mutate(Mouse_ID = factor(Mouse_ID)) %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS FOS+GOS"))) %>% 
  group_by(sex) %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_pups_FOSGOS_sex_w5
write_xlsx(tukey_pups_FOSGOS_sex_w5, "tukey_pups_FOSGOS_sex_w5.xlsx")

##Tukeys pups FOS+GOS grouped by sex w10 only
tukey_pups_FOSGOS_sex_w10 = metadata_pup %>% 
  filter(diet != "HFHS B. longum APC1472") %>% 
  filter(timepoint != "w5") %>%
  #mutate(Mouse_ID = factor(Mouse_ID)) %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS FOS+GOS"))) %>% 
  group_by(sex) %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_pups_FOSGOS_sex_w10
write_xlsx(tukey_pups_FOSGOS_sex_w10, "tukey_pups_FOSGOS_sex_w10.xlsx")

##Tukeys pups FOS+GOS grouped by timepoint
tukey_pups_FOSGOS_timepoint = metadata_pup %>% 
  filter(diet != "HFHS B. longum APC1472") %>% 
  #mutate(Mouse_ID = factor(Mouse_ID)) %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS FOS+GOS"))) %>% 
  group_by(timepoint) %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_pups_FOSGOS_timepoint
write_xlsx(tukey_pups_FOSGOS_timepoint, "tukey_pups_FOSGOS_timepoint.xlsx")

#Tukeys pups 1472 grouped by sex
tukey_pups_1472_sex = metadata_pup %>% 
  filter(diet != "HFHS FOS+GOS") %>% 
  #mutate(Mouse_ID = factor(Mouse_ID)) %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS B. longum APC1472"))) %>% 
  group_by(sex) %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_pups_1472_sex
write_xlsx(tukey_pups_1472_sex, "tukey_pups_1472_sex.xlsx")

#Tukeys pups 1472 grouped by sex w5
tukey_pups_1472_sex_w5 = metadata_pup %>% 
  filter(diet != "HFHS FOS+GOS") %>% 
  filter(timepoint != "w10") %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS B. longum APC1472"))) %>% 
  group_by(sex) %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_pups_1472_sex_w5
write_xlsx(tukey_pups_1472_sex_w5, "tukey_pups_1472_sex_w5.xlsx")

#Tukeys pups 1472 grouped by sex w10
tukey_pups_1472_sex_w10 = metadata_pup %>% 
  filter(diet != "HFHS FOS+GOS") %>% 
  filter(timepoint != "w5") %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS B. longum APC1472"))) %>% 
  group_by(sex) %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_pups_1472_sex_w10
write_xlsx(tukey_pups_1472_sex_w10, "tukey_pups_1472_sex_w10.xlsx")

#Tukeys pups 1472 grouped by timepoint
tukey_pups_1472_timepoint = metadata_pup %>% 
  filter(diet != "HFHS FOS+GOS") %>% 
  #mutate(Mouse_ID = factor(Mouse_ID)) %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS B. longum APC1472"))) %>% 
  group_by(timepoint) %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_pups_1472_timepoint
write_xlsx(tukey_pups_1472_timepoint, "tukey_pups_1472_timepoint.xlsx")



##########dams
# Extract sample data and abundance data
metadata_dam <- data.frame(sample_data(PSr_dam))  # Extract sample metadata
otu_data_dam <- data.frame(otu_table(PSr_dam))

rownames(metadata_dam) <- sample_names(PSr_dam)

#####alpha diversity
div.all.1_dam <- estimate_richness(ps_dam, measures = "Shannon")
metadata_dam$shan <- div.all.1_dam

rownames(metadata_dam) == rownames(div.all.1_dam)

##Summary lm-FOS+GOS
mat_restored_FOSGOS  = cbind(c(1/2, -1, 1/2))
rownames(mat_restored_FOSGOS) <- c("Control Vehicle","HFHS Vehicle","HFHS FOS+GOS")
colnames(mat_restored_FOSGOS) <- c("Restored")

lm_summary_dams_FOSGOS <- metadata_dam %>%
  filter(diet != "HFHS B. longum APC1472") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS"))) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_FOSGOS),
       data = . 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_summary_dams_FOSGOS
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Stats")
write_xlsx(lm_summary_dams_FOSGOS, "lm_summary_dams_FOSGOS.xlsx")

##Summary lm-1472
mat_restored_1472  = cbind(c(1/2, -1, 1/2))
rownames(mat_restored_1472) <- c("Control Vehicle","HFHS Vehicle","HFHS B. longum APC1472")
colnames(mat_restored_1472) <- c("Restored")

lm_summary_dams_1472 <- metadata_dam %>%
  filter(diet != "HFHS FOS+GOS") %>%
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS B. longum APC1472"))) %>%
  reframe(
    lm(shan$Shannon ~  diet, 
       contrasts = list(diet = mat_restored_1472),
       data = . 
    ) %>% 
      car::Anova(., type = 3) %>% 
      tidy()) %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p.adjusted = p.adjust(p = p.value, "BH")) %>%
  ungroup
lm_summary_dams_1472
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Stats")
write_xlsx(lm_summary_dams_1472, "lm_summary_dams_1472.xlsx")

###Tukey dams 1472
tukey_dams_1472 = metadata_dam %>% 
  filter(diet != "HFHS FOS+GOS") %>% 
  #mutate(Mouse_ID = factor(Mouse_ID)) %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS B. longum APC1472"))) %>% 
  group_by() %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_dams_1472
write_xlsx(tukey_dams_1472, "tukey_dams_1472.xlsx")

##Tukey FOS GOS
tukey_dams_FOSGOS = metadata_dam %>% 
  filter(diet != "HFHS B. longum APC1472") %>% 
  #mutate(Mouse_ID = factor(Mouse_ID)) %>% 
  mutate(diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle","HFHS FOS+GOS"))) %>% 
  group_by() %>% 
  reframe(
    lm(shan$Shannon ~ diet, 
       data = . 
    ) %>% aov() %>% TukeyHSD() %>% tidy()
  ) %>% 
  ungroup() %>%
  filter(term != "(Intercept)")
tukey_dams_FOSGOS
write_xlsx(tukey_dams_FOSGOS, "tukey_dams_FOSGOS.xlsx")
