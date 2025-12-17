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
library(writexl)
library(janitor)
library(tibble)
library(purrr)
library(rstatix)
library(data.table)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(decontam)

## Make phyloseq object
setwd("C:/Users/DHedayatpour/OneDrive - University College Cork/ADATA onedrive/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items")
getwd()

samdf <- read_excel("Metadata_BINC_maxee24_ps.xlsx")

sam2 <- samdf %>% remove_rownames %>% column_to_rownames(var="sample")
sam2$sample <- rownames(sam2)

# split my diet column into diet and treatment
sam2 <- sam2 %>%
  separate(diet, into = c("diet", "treatment"), sep = " ", extra = "merge") %>%
  mutate(diet = factor(diet, levels = c("Control", "HFHS")),
         treatment = factor(treatment, levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         sex = factor(sex, levels = c("F", "M")),
         timepoint = factor(timepoint, levels = c("w5", "w10")))

sam2$sex <- factor(sam2$sex, levels = c("M", "F"), labels = c("Male", "Female"))
sam2$timepoint <- factor(sam2$timepoint, levels = c("w5", "w10"), labels = c("week 5", "week 10"))

##load the tax table
taxa <- read.csv("taxa_maxee24_ps.csv", sep = ",", row.names = 1)
taxa_mat <- as.matrix(taxa)

##load the otu table
seqtab.nochim <- read.csv("seqtab_nochim_maxee24_ps.csv", sep = ",", row.names = 1)
seqtab.nochim_t <- t(seqtab.nochim)

seq_mat <- as.matrix(seqtab.nochim)

###make a phyloseq object 

ps <- phyloseq(otu_table(seq_mat, taxa_are_rows=FALSE), 
               sample_data(sam2), 
               tax_table(taxa_mat))

ps

### decontam
#Put sample_data into a ggplot-friendly: data.frame
df <- data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
#note the new column added: LibrarySize
#now, order df by LibrarySize, small to big
df <- df[order(df$LibrarySize),]
df
df$Index <- seq(nrow(df))
df
tail(df)
ggplot(data=df, aes(x=Index, y=LibrarySize)) +
  geom_point()#color=Sample_or_control)) + geom_point()

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_control == "control"
contamdf.prev2 <- isContaminant(ps, method="prevalence", neg="is.neg")
#contamdf.prev2 <- isContaminant(ps, method="frequency", conc="quant_reading")
table(contamdf.prev2$contaminant)
# FALSE 2716
head(which(contamdf.prev2$contaminant))
# integer(0)
##now adding threshold = 0.5
contamdf.prev2 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev2$contaminant)
# FALSE 2716

ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_control == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$is.neg == "FALSE", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                     contaminant=contamdf.prev2$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

pps <- prune_taxa(!contamdf.prev2$contaminant, ps)  

pps # no contaminants

##remove control too
look_at_this <- psmelt(ps)
look_at_this_neg <- look_at_this %>%
  filter(Sample == "negblank") # I have 0 abundance of any OTU in my negblank.
pps_nocontam <- subset_samples(pps, Sample_or_control != "control")

##remove non-bacteria also
pps_bacteria <- pps_nocontam %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast" &
      Phylum != "Cyanobacteria/Chloroplast"
  )

pps_bacteria

##separate the pups and dams, and make separate ps based on age column
ps_pup <- subset_samples(pps_bacteria, sample_data(pps_bacteria)$age =="pup")
ps_dam <- subset_samples(pps_bacteria, sample_data(pps_bacteria)$age =="dam")

##Make relative abundances
PSr_pup <- transform_sample_counts(ps_pup, function(x) x/sum(x))
PSr_dam <- transform_sample_counts(ps_dam, function(x) x/sum(x))

############## Bifidobacterium species
############## Genus
############## Family
############## Order
############## Phylum
############## Alpha diversity shannon
############## Alpha diversity simpson
############## Alpha diversity chao1
############## Beta diversity
############## Bacillota:Bacteroidota ratio

####### now to do bif

##Filter for Bif species relative
bif_relative_pup <- PSr_pup %>%
  tax_glom(taxrank = "Species") %>%
  psmelt() %>%
  filter(Genus == "Bifidobacterium") %>%
  mutate(Genus_species = paste(Genus, Species, sep = " "))

##Filter for bif absolute
bif_absolute_pup <- ps_pup %>%
  tax_glom(taxrank = "Species") %>%
  psmelt() %>%
  filter(Genus == "Bifidobacterium") %>%
  mutate(Genus_species = paste(Genus, Species, sep = " "))

##Part 1 asks: how does the interaction between diet and treatment affect the microbiome. add cage and cohort as random effects.
## random effect like dam_cage means the linear model now factors that pups born from the same cage (dam) may influence the result.
bif_species_pups_relabund_diet_treatment_interaction <- bif_relative_pup %>%
  select(Genus_species, cohort, sex, dam_cage, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
   
    #pull(treatment) %>% unique() %>%
    #pull(diet) %>% unique() %>%
   
  mutate(treatment = factor(x = treatment,
                              levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                         levels = c("Control", "HFHS")),
         cohort = factor(cohort),
         dam_cage = factor(dam_cage)
           ) %>%
   
  group_by(genus_species) %>%
  nest() %>%
  mutate(lm = purrr::map(.x = data,
                           .f = ~ if(var(.x$abundance) == 0) {
                             tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
                           } else {
                             lmerTest::lmer(abundance ~ diet * treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
                               broom.mixed::tidy()} 
                           )) %>%
  select(!data) %>%
  unnest(lm) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate))

write.csv(bif_species_pups_relabund_diet_treatment_interaction, "bif_species_pups_relabund_diet_treatment_interaction.csv", row.names = FALSE)


###Diet---
## Part 2: tests if diet (control or HFHS) has an effect in vehicle treatment group (no APC1472 or FOS+GOS) per bif species, sex, and timepoint.
diet_bif_species_relabund_pups <- bif_relative_pup %>%
  select(Genus_species, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                              levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                         levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
        cohort = factor(cohort), 
        sex = factor(sex, levels = c("Female", "Male")),
        timepoint = factor(timepoint, levels = c("week 5", "week 10"))
        ) %>%
  filter(treatment == "Vehicle") %>%
  group_by(genus_species, sex, timepoint) %>%
  nest() %>%
  mutate(diet_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$abundance) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(abundance ~ diet + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>% 
  unnest(diet_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "diet")
 
###Treatment
## Part 3: tests if treatment (vehicle, APC1472, or FOS+GOS) has an effect on HFHS diet group, per bif species, for sex and timepoint.
  treatment_bif_species_relabund_pups <- bif_relative_pup %>%
    select(Genus_species, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    mutate(treatment = factor(x = treatment,
                              levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
           diet = factor(x = diet,
                         levels = c("Control", "HFHS")),
           dam_cage = factor(dam_cage),
           cohort = factor(cohort), 
           sex = factor(sex, levels = c("Female", "Male")),
           timepoint = factor(timepoint, levels = c("week 5", "week 10"))
    ) %>%
    filter(diet == "HFHS") %>%
    group_by(genus_species, sex, timepoint) %>%
    nest() %>%
    mutate(treatment_lmer = purrr::map(.x = data, .f = ~ {
      if(var(.x$abundance) == 0) {
        tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
      } else {
        lmerTest::lmer(abundance ~ treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
          broom.mixed::tidy()
      }
    })
    ) %>%
    select(-data) %>%  
    unnest(treatment_lmer) %>%
    mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
    filter(term != "(Intercept)" & !is.na(estimate)) %>%
    mutate(comparison_type = "treatment")
 
 ## Part 4: combining results from part 2 and 3
combined_diet_treatment_bif_species_relabund_pups <- bind_rows(diet_bif_species_relabund_pups, treatment_bif_species_relabund_pups)
 
write.csv(combined_diet_treatment_bif_species_relabund_pups,
            "combined_diet_treatment_bif_species_relabund_pups.csv",
            row.names = FALSE) 

 ###Restoration
 ## Part 5: Checks whether treatment actually reversed the effects of HFHS i.e restore microbiome back to control
  diet_treatment_bif_species_relabund_pups_restoration <- combined_diet_treatment_bif_species_relabund_pups %>%
    filter(
      (comparison_type == "diet" & term == "dietHFHS") | # Benjamin's script has group1/2, whereas because I did lmer to account for random variable, I have term instead of group1/2.
        (comparison_type == "treatment" & term %in% c("treatmentFOS+GOS", "treatmentB. longum APC1472"))
    ) %>%
    group_by(genus_species, sex, timepoint) %>%
    summarise(
      statistic_diet = first(statistic[term == "dietHFHS"], default = NA),
      p.adj_diet = first(p.adj[term == "dietHFHS"], default = NA), 
      statistic_treatment_FOSGOS = first(statistic[term == "treatmentFOS+GOS"], default = NA),
      p.adj_treatment_FOSGOS = first(p.adj[term == "treatmentFOS+GOS"], default = NA),
      statistic_treatment_1472 = first(statistic[term == "treatmentB. longum APC1472"], default = NA),
      p.adj_treatment_1472 = first(p.adj[term == "treatmentB. longum APC1472"], default = NA)
    ) %>%
    mutate(
      restoration_FOSGOS = if_else(
        !is.na(p.adj_diet) & !is.na(p.adj_treatment_FOSGOS) &  
          p.adj_diet < 0.05 & p.adj_treatment_FOSGOS < 0.05 &  
          sign(statistic_diet) != sign(statistic_treatment_FOSGOS),  
        "restored", "not restored"
      ),
      restoration_1472 = if_else(
        !is.na(p.adj_diet) & !is.na(p.adj_treatment_1472) &  
          p.adj_diet < 0.05 & p.adj_treatment_1472 < 0.05 &  
          sign(statistic_diet) != sign(statistic_treatment_1472),  
        "restored", "not restored"
      )
    )
  write.csv(diet_treatment_bif_species_relabund_pups_restoration, "diet_treatment_bif_species_relabund_pups_restoration.csv", row.names = FALSE)

  
  ## why is longum only in male week 10:
  # for longum in Fw5, Fw10, and Mw10, the relative abundance was so low that there was no variation available for modelling.
  
######## Genus
##Filter for genus relative
genus_relative_pup <- PSr_pup %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()

##Part 1 asks: how does the interaction between diet and treatment affect the microbiome. 
genus_pups_relabund_diet_treatment_interaction <- genus_relative_pup %>%
  select(Genus, cohort, sex, dam_cage, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  
  #pull(treatment) %>% unique() %>%
  #pull(diet) %>% unique() %>%
  
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         cohort = factor(cohort),
         dam_cage = factor(dam_cage)
  ) %>%
  
  group_by(genus) %>%
  nest() %>%
  mutate(lm = purrr::map(.x = data,
                         .f = ~ if(var(.x$abundance) == 0) {
                           tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
                         } else {
                           lmerTest::lmer(abundance ~ diet * treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
                             broom.mixed::tidy()} 
  )) %>%
  select(!data) %>%
  unnest(lm) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate))

write.csv(genus_pups_relabund_diet_treatment_interaction, "genus_pups_relabund_diet_treatment_interaction.csv", row.names = FALSE)

###Diet---
## Part 2: tests if diet (control or HFHS) has an effect in vehicle treatment group (no APC1472 or FOS+GOS) per genus, sex, and timepoint.
diet_genus_relabund_pups <- genus_relative_pup %>%
  select(Genus, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(treatment == "Vehicle") %>%
  group_by(genus, sex, timepoint) %>%
  nest() %>%
  mutate(diet_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$abundance) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(abundance ~ diet + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(diet_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "diet")

###Treatment
## Part 3: tests if treatment (vehicle, APC1472, or FOS+GOS) has an effect on HFHS diet group, per genus, for sex and timepoint.
treatment_genus_relabund_pups <- genus_relative_pup %>%
  select(Genus, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(diet == "HFHS") %>%
  group_by(genus, sex, timepoint) %>%
  nest() %>%
  mutate(treatment_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$abundance) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(abundance ~ treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(treatment_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "treatment")

## Part 4: combining results from part 2 and 3
combined_diet_treatment_genus_relabund_pups <- bind_rows(diet_genus_relabund_pups, treatment_genus_relabund_pups)

write.csv(combined_diet_treatment_genus_relabund_pups,
          "combined_diet_treatment_genus_relabund_pups.csv",
          row.names = FALSE) 

###Restoration
## Part 5: Checks whether treatment actually reversed the effects of HFHS i.e restore microbiome back to control
diet_treatment_genus_relabund_pups_restoration <- combined_diet_treatment_genus_relabund_pups %>%
  filter(
    (comparison_type == "diet" & term == "dietHFHS") | 
      (comparison_type == "treatment" & term %in% c("treatmentFOS+GOS", "treatmentB. longum APC1472"))
  ) %>%
  group_by(genus, sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[term == "dietHFHS"], default = NA),
    p.adj_diet = first(p.adj[term == "dietHFHS"], default = NA), 
    statistic_treatment_FOSGOS = first(statistic[term == "treatmentFOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[term == "treatmentFOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[term == "treatmentB. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[term == "treatmentB. longum APC1472"], default = NA)
  ) %>%
  mutate(
    restoration_FOSGOS = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_FOSGOS) &  
        p.adj_diet < 0.05 & p.adj_treatment_FOSGOS < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_FOSGOS),  
      "restored", "not restored"
    ),
    restoration_1472 = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_1472) &  
        p.adj_diet < 0.05 & p.adj_treatment_1472 < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_1472),  
      "restored", "not restored"
    )
  )
write.csv(diet_treatment_genus_relabund_pups_restoration, "diet_treatment_genus_relabund_pups_restoration.csv", row.names = FALSE)

######## Family
##Filter for family relative
family_relative_pup <- PSr_pup %>%
  tax_glom(taxrank = "Family") %>%
  psmelt()

##Part 1 asks: how does the interaction between diet and treatment affect the microbiome. 
family_pups_relabund_diet_treatment_interaction <- family_relative_pup %>%
  select(Family, cohort, sex, dam_cage, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  
  #pull(treatment) %>% unique() %>%
  #pull(diet) %>% unique() %>%
  
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         cohort = factor(cohort),
         dam_cage = factor(dam_cage)
  ) %>%
  
  group_by(family) %>%
  nest() %>%
  mutate(lm = purrr::map(.x = data,
                         .f = ~ if(var(.x$abundance) == 0) {
                           tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
                         } else {
                           lmerTest::lmer(abundance ~ diet * treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
                             broom.mixed::tidy()} 
  )) %>%
  select(!data) %>%
  unnest(lm) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate))

write.csv(family_pups_relabund_diet_treatment_interaction, "family_pups_relabund_diet_treatment_interaction.csv", row.names = FALSE)

###Diet---
## Part 2: tests if diet (control or HFHS) has an effect in vehicle treatment group (no APC1472 or FOS+GOS) per family, sex, and timepoint.
diet_family_relabund_pups <- family_relative_pup %>%
  select(Family, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(treatment == "Vehicle") %>%
  group_by(family, sex, timepoint) %>%
  nest() %>%
  mutate(diet_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$abundance) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(abundance ~ diet + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(diet_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "diet")

###Treatment
## Part 3: tests if treatment (vehicle, APC1472, or FOS+GOS) has an effect on HFHS diet group, per family, for sex and timepoint.
treatment_family_relabund_pups <- family_relative_pup %>%
  select(Family, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(diet == "HFHS") %>%
  group_by(family, sex, timepoint) %>%
  nest() %>%
  mutate(treatment_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$abundance) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(abundance ~ treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(treatment_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "treatment")

## Part 4: combining results from part 2 and 3
combined_diet_treatment_family_relabund_pups <- bind_rows(diet_family_relabund_pups, treatment_family_relabund_pups)

write.csv(combined_diet_treatment_family_relabund_pups,
          "combined_diet_treatment_family_relabund_pups.csv",
          row.names = FALSE) 

###Restoration
## Part 5: Checks whether treatment actually reversed the effects of HFHS i.e restore microbiome back to control
diet_treatment_family_relabund_pups_restoration <- combined_diet_treatment_family_relabund_pups %>%
  filter(
    (comparison_type == "diet" & term == "dietHFHS") | 
      (comparison_type == "treatment" & term %in% c("treatmentFOS+GOS", "treatmentB. longum APC1472"))
  ) %>%
  group_by(family, sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[term == "dietHFHS"], default = NA),
    p.adj_diet = first(p.adj[term == "dietHFHS"], default = NA), 
    statistic_treatment_FOSGOS = first(statistic[term == "treatmentFOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[term == "treatmentFOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[term == "treatmentB. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[term == "treatmentB. longum APC1472"], default = NA)
  ) %>%
  mutate(
    restoration_FOSGOS = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_FOSGOS) &  
        p.adj_diet < 0.05 & p.adj_treatment_FOSGOS < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_FOSGOS),  
      "restored", "not restored"
    ),
    restoration_1472 = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_1472) &  
        p.adj_diet < 0.05 & p.adj_treatment_1472 < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_1472),  
      "restored", "not restored"
    )
  )
write.csv(diet_treatment_family_relabund_pups_restoration, "diet_treatment_family_relabund_pups_restoration.csv", row.names = FALSE)

######## Order
##Filter for order relative
order_relative_pup <- PSr_pup %>%
  tax_glom(taxrank = "Order") %>%
  psmelt()

##Part 1 asks: how does the interaction between diet and treatment affect the microbiome.
order_pups_relabund_diet_treatment_interaction <- order_relative_pup %>%
  select(Order, cohort, sex, dam_cage, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         cohort = factor(cohort),
         dam_cage = factor(dam_cage)
  ) %>%
  
  group_by(order) %>%
  nest() %>%
  mutate(lm = purrr::map(.x = data,
                         .f = ~ if(var(.x$abundance) == 0) {
                           tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
                         } else {
                           lmerTest::lmer(abundance ~ diet * treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
                             broom.mixed::tidy()} 
  )) %>%
  select(!data) %>%
  unnest(lm) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate))

write.csv(order_pups_relabund_diet_treatment_interaction, "order_pups_relabund_diet_treatment_interaction.csv", row.names = FALSE)

###Diet---
## Part 2: tests if diet (control or HFHS) has an effect in vehicle treatment group (no APC1472 or FOS+GOS) per order, sex, and timepoint.
diet_order_relabund_pups <- order_relative_pup %>%
  select(Order, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(treatment == "Vehicle") %>%
  group_by(order, sex, timepoint) %>%
  nest() %>%
  mutate(diet_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$abundance) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(abundance ~ diet + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(diet_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "diet")

###Treatment
## Part 3: tests if treatment (vehicle, APC1472, or FOS+GOS) has an effect on HFHS diet group, per order, for sex and timepoint.
treatment_order_relabund_pups <- order_relative_pup %>%
  select(Order, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(diet == "HFHS") %>%
  group_by(order, sex, timepoint) %>%
  nest() %>%
  mutate(treatment_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$abundance) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(abundance ~ treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>% 
  unnest(treatment_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "treatment")

## Part 4: combining results from part 2 and 3
combined_diet_treatment_order_relabund_pups <- bind_rows(diet_order_relabund_pups, treatment_order_relabund_pups)

write.csv(combined_diet_treatment_order_relabund_pups,
          "combined_diet_treatment_order_relabund_pups.csv",
          row.names = FALSE) 

###Restoration
## Part 5: Checks whether treatment actually reversed the effects of HFHS i.e restore microbiome back to control
diet_treatment_order_relabund_pups_restoration <- combined_diet_treatment_order_relabund_pups %>%
  filter(
    (comparison_type == "diet" & term == "dietHFHS") | 
      (comparison_type == "treatment" & term %in% c("treatmentFOS+GOS", "treatmentB. longum APC1472"))
  ) %>%
  group_by(order, sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[term == "dietHFHS"], default = NA),
    p.adj_diet = first(p.adj[term == "dietHFHS"], default = NA), 
    statistic_treatment_FOSGOS = first(statistic[term == "treatmentFOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[term == "treatmentFOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[term == "treatmentB. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[term == "treatmentB. longum APC1472"], default = NA)
  ) %>%
  mutate(
    restoration_FOSGOS = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_FOSGOS) &  
        p.adj_diet < 0.05 & p.adj_treatment_FOSGOS < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_FOSGOS),  
      "restored", "not restored"
    ),
    restoration_1472 = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_1472) &  
        p.adj_diet < 0.05 & p.adj_treatment_1472 < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_1472),  
      "restored", "not restored"
    )
  )
write.csv(diet_treatment_order_relabund_pups_restoration, "diet_treatment_order_relabund_pups_restoration.csv", row.names = FALSE)

######## Phylum
##Filter for phylum relative
phylum_relative_pup <- PSr_pup %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt()

##Part 1 asks: how does the interaction between diet and treatment affect the microbiome.
phylum_pups_relabund_diet_treatment_interaction <- phylum_relative_pup %>%
  select(Phylum, cohort, sex, dam_cage, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         cohort = factor(cohort),
         dam_cage = factor(dam_cage)
  ) %>%
  
  group_by(phylum) %>%
  nest() %>%
  mutate(lm = purrr::map(.x = data,
                         .f = ~ if(var(.x$abundance) == 0) {
                           tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
                         } else {
                           lmerTest::lmer(abundance ~ diet * treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
                             broom.mixed::tidy()} 
  )) %>%
  select(!data) %>%
  unnest(lm) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate))

write.csv(phylum_pups_relabund_diet_treatment_interaction, "phylum_pups_relabund_diet_treatment_interaction.csv", row.names = FALSE)

###Diet---
## Part 2: tests if diet (control or HFHS) has an effect in vehicle treatment group (no APC1472 or FOS+GOS) per phylum, sex, and timepoint.
diet_phylum_relabund_pups <- phylum_relative_pup %>%
  select(Phylum, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(treatment == "Vehicle") %>%
  group_by(phylum, sex, timepoint) %>%
  nest() %>%
  mutate(diet_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$abundance) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(abundance ~ diet + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(diet_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "diet")

###Treatment
## Part 3: tests if treatment (vehicle, APC1472, or FOS+GOS) has an effect on HFHS diet group, per phylum, for sex and timepoint.
treatment_phylum_relabund_pups <- phylum_relative_pup %>%
  select(Phylum, cohort, dam_cage, sex, diet, timepoint, treatment, Sample, Abundance) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(diet == "HFHS") %>%
  group_by(phylum, sex, timepoint) %>%
  nest() %>%
  mutate(treatment_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$abundance) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(abundance ~ treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>% 
  unnest(treatment_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "treatment")

## Part 4: combining results from part 2 and 3
combined_diet_treatment_phylum_relabund_pups <- bind_rows(diet_phylum_relabund_pups, treatment_phylum_relabund_pups)

write.csv(combined_diet_treatment_phylum_relabund_pups,
          "combined_diet_treatment_phylum_relabund_pups.csv",
          row.names = FALSE) 

###Restoration
## Part 5: Checks whether treatment actually reversed the effects of HFHS i.e restore microbiome back to control
diet_treatment_phylum_relabund_pups_restoration <- combined_diet_treatment_phylum_relabund_pups %>%
  filter(
    (comparison_type == "diet" & term == "dietHFHS") | 
      (comparison_type == "treatment" & term %in% c("treatmentFOS+GOS", "treatmentB. longum APC1472"))
  ) %>%
  group_by(phylum, sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[term == "dietHFHS"], default = NA),
    p.adj_diet = first(p.adj[term == "dietHFHS"], default = NA), 
    statistic_treatment_FOSGOS = first(statistic[term == "treatmentFOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[term == "treatmentFOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[term == "treatmentB. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[term == "treatmentB. longum APC1472"], default = NA)
  ) %>%
  mutate(
    restoration_FOSGOS = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_FOSGOS) &  
        p.adj_diet < 0.05 & p.adj_treatment_FOSGOS < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_FOSGOS),  
      "restored", "not restored"
    ),
    restoration_1472 = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_1472) &  
        p.adj_diet < 0.05 & p.adj_treatment_1472 < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_1472),  
      "restored", "not restored"
    )
  )
write.csv(diet_treatment_phylum_relabund_pups_restoration, "diet_treatment_phylum_relabund_pups_restoration.csv", row.names = FALSE)

###### Beta diversity
##Calculate beta diversity distance. Use relative abundances.
DistBC <- distance(PSr_pup, method = "bray")

##Perform PCoA on each distance matrix
ordBC <- ordinate(PSr_pup, method = "PCoA", distance = DistBC)
ordBC$values$Eigenvalues
eig_values <- ordBC$values$Eigenvalues
percent_explained <- (eig_values / sum(eig_values)) * 100
percent_explained[1:5] # axis 1 accounts for 20-40 percent variation
plot_scree(ordBC, "Scree Plot: Bray-Curtis PSr_pup")  

#convert phyloseq sample data dataframe to a data table
sdt <- data.table(as(sample_data(PSr_pup), "data.frame"),
                  keep.rownames = TRUE) # rn (rownames) matches the 196 sample names
setnames(sdt, "rn", "SampleID") # rn is now called SampleID

##Join sample data and ordination axes together in one data.table
ordBCdt <- data.table(as.data.frame(ordBC$vectors), keep.rownames = TRUE)
setnames(ordBCdt, "rn", "SampleID")
setkey(ordBCdt, SampleID)
setkey(sdt, SampleID)
#ordBCsdt <- ordBCdt[sdt]
ordBCsdt <- sdt[ordBCdt]
#setorder(ordBCsdt, timepoint)
str(ordBCsdt)

###Part 1: how does the interaction between diet and treatment affect the beta diversity With cage and cohort as random effects.
beta_pups_relabund_diet_treatment_interaction <- ordBCsdt %>%
  select(cohort, sex, dam_cage, diet, timepoint, treatment, sample, Axis.1, Axis.2, Axis.3) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10")),
         cohort = factor(cohort),
         dam_cage = factor(dam_cage)
  ) %>%
  
  #group_by(axis_1) %>%
  nest() %>%
  mutate(lm = purrr::map(.x = data,
                         .f = ~ if (var(.x$axis_1) == 0) {
                           tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
                         } else {
                           lmerTest::lmer(axis_1 ~ diet * treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
                             broom.mixed::tidy()
                         }
  )) %>%
  select(!data) %>%
  unnest(lm) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate))

write.csv(beta_pups_relabund_diet_treatment_interaction, "beta_pups_relabund_diet_treatment_interaction.csv", row.names = FALSE)

###Diet---
## Part 2: tests if diet (control or HFHS) has an effect in vehicle treatment group (no APC1472 or FOS+GOS) beta diversity, sex, and timepoint.
diet_beta_relabund_pups <- ordBCsdt %>%
  select(cohort, dam_cage, sex, diet, timepoint, treatment, sample, Axis.1) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(treatment == "Vehicle") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(diet_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$axis_1) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(axis_1 ~ diet + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(diet_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "diet")

###Treatment
## Part 3: tests if treatment (vehicle, APC1472, or FOS+GOS) has an effect on HFHS diet group, per beta diversity, for sex and timepoint.
treatment_beta_relabund_pups <- ordBCsdt %>%
  select(cohort, dam_cage, sex, diet, timepoint, treatment, sample, Axis.1) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(diet == "HFHS") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(treatment_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$axis_1) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(axis_1 ~ treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(treatment_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "treatment")

## Part 4: combining results from part 2 and 3
combined_diet_treatment_beta_relabund_pups <- bind_rows(diet_beta_relabund_pups, treatment_beta_relabund_pups)

write.csv(combined_diet_treatment_beta_relabund_pups,
          "combined_diet_treatment_beta_relabund_pups.csv",
          row.names = FALSE)

###Restoration
## Part 5: Checks whether treatment actually reversed the effects of HFHS i.e restore microbiome back to control
diet_treatment_beta_relabund_pups_restoration <- combined_diet_treatment_beta_relabund_pups %>%
  filter(
    (comparison_type == "diet" & term == "dietHFHS") | 
      (comparison_type == "treatment" & term %in% c("treatmentFOS+GOS", "treatmentB. longum APC1472"))
  ) %>%
  group_by(sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[term == "dietHFHS"], default = NA),
    p.adj_diet = first(p.adj[term == "dietHFHS"], default = NA), 
    statistic_treatment_FOSGOS = first(statistic[term == "treatmentFOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[term == "treatmentFOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[term == "treatmentB. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[term == "treatmentB. longum APC1472"], default = NA)
  ) %>%
  mutate(
    restoration_FOSGOS = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_FOSGOS) &  
        p.adj_diet < 0.05 & p.adj_treatment_FOSGOS < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_FOSGOS),  
      "restored", "not restored"
    ),
    restoration_1472 = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_1472) &  
        p.adj_diet < 0.05 & p.adj_treatment_1472 < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_1472),  
      "restored", "not restored"
    )
  )
write.csv(diet_treatment_beta_relabund_pups_restoration, "diet_treatment_beta_relabund_pups_restoration.csv", row.names = FALSE)

############# Alpha Diversity
alpha_diversity <- estimate_richness(ps_pup, measures = c("Shannon", "Chao1", "Simpson")) 
rownames(alpha_diversity) <- gsub("^X", "", rownames(alpha_diversity)) # remove pesky X

# Merge with sample metadata VERY IMPORTANT STEP
meta <- sample_data(ps_pup)
alpha_diversity <- merge(meta, alpha_diversity, by = "row.names")
row.names(meta)
row.names(alpha_diversity)

###Part 1: how does the interaction between diet and treatment affect the alpha diversity With cage and cohort as random effects.
### Shannon
alpha_pups_diet_treatment_interaction <- alpha_diversity %>%
  select(cohort, sex, dam_cage, diet, timepoint, treatment, sample, Shannon, Simpson, Chao1) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  pivot_longer(cols = c(shannon, simpson, chao1), 
               names_to = "metric", values_to = "value") %>% # allows for the three metrics to fit the models at once
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10")),
         cohort = factor(cohort),
         dam_cage = factor(dam_cage)
  ) %>%
  group_by(metric) %>%
  nest() %>%
  mutate(lm = purrr::map(.x = data,
                         .f = ~ if (var(.x$value) == 0) {
                           tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
                         } else {
                           lmerTest::lmer(value ~ diet * treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
                             broom.mixed::tidy()
                         }
  )) %>%
  select(!data) %>%
  unnest(lm) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate))

write.csv(alpha_pups_diet_treatment_interaction, "alpha_pups_diet_treatment_interaction.csv", row.names = FALSE)

###Diet---
## Part 2: tests if diet (control or HFHS) has an effect in vehicle treatment group (no APC1472 or FOS+GOS) per alpha diversity, sex, and timepoint.
diet_alpha_pups <- alpha_diversity %>%
  select(cohort, dam_cage, sex, diet, timepoint, treatment, sample, Shannon, Simpson, Chao1) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  pivot_longer(cols = c(shannon, simpson, chao1),
               names_to = "metric", values_to = "value") %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(treatment == "Vehicle") %>%
  group_by(metric, sex, timepoint) %>%
  nest() %>%
  mutate(diet_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$value) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(value ~ diet + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(diet_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "diet")

###Treatment
## Part 3: tests if treatment (vehicle, APC1472, or FOS+GOS) has an effect on HFHS diet group, per alpha diversity, for sex and timepoint.
treatment_alpha_pups <- alpha_diversity %>%
  select(cohort, dam_cage, sex, diet, timepoint, treatment, sample, Shannon, Simpson, Chao1) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  pivot_longer(cols = c(shannon, simpson, chao1),
               names_to = "metric", values_to = "value") %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(diet == "HFHS") %>%
  group_by(metric, sex, timepoint) %>%
  nest() %>%
  mutate(treatment_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$value) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(value ~ treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>% 
  unnest(treatment_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "treatment")

## Part 4: combining results from part 2 and 3
combined_diet_treatment_alpha_pups <- bind_rows(diet_alpha_pups, treatment_alpha_pups)

write.csv(combined_diet_treatment_alpha_pups,
          "combined_diet_treatment_alpha_pups.csv",
          row.names = FALSE)

###Restoration
## Part 5: Checks whether treatment actually reversed the effects of HFHS i.e restore microbiome back to control
diet_treatment_alpha_pups_restoration <- combined_diet_treatment_alpha_pups %>%
  filter(
    (comparison_type == "diet" & term == "dietHFHS") | 
      (comparison_type == "treatment" & term %in% c("treatmentFOS+GOS", "treatmentB. longum APC1472"))
  ) %>%
  group_by(metric, sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[term == "dietHFHS"], default = NA),
    p.adj_diet = first(p.adj[term == "dietHFHS"], default = NA), 
    statistic_treatment_FOSGOS = first(statistic[term == "treatmentFOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[term == "treatmentFOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[term == "treatmentB. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[term == "treatmentB. longum APC1472"], default = NA)
  ) %>%
  mutate(
    restoration_FOSGOS = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_FOSGOS) &  
        p.adj_diet < 0.05 & p.adj_treatment_FOSGOS < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_FOSGOS),  
      "restored", "not restored"
    ),
    restoration_1472 = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_1472) &  
        p.adj_diet < 0.05 & p.adj_treatment_1472 < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_1472),  
      "restored", "not restored"
    )
  )
write.csv(diet_treatment_alpha_pups_restoration, "diet_treatment_alpha_pups_restoration.csv", row.names = FALSE)

############# Bacillota:Bacteroidota ratio
ps_pup_phylum <- PSr_pup %>%
  psmelt() %>%
  mutate(Phylum = ifelse(Phylum == "Firmicutes", "Bacillota", Phylum)) # rename firmicutes

bac_pup_ratio <- ps_pup_phylum %>%
  filter(Phylum %in% c("Bacillota", "Bacteroidota")) %>%
  group_by(Sample, diet, timepoint, sex, treatment, cohort, dam_cage, Phylum) %>%
  summarise(TotalRelAbundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Phylum, values_from = TotalRelAbundance, values_fill = 0) %>%
  mutate(Bac_Bact_ratio = Bacillota / Bacteroidota) # Bac_Bact_ratio as variable

###Part 1: how does the interaction between diet and treatment affect the ratio With cage and cohort as random effects.
ratio_pups_diet_treatment_interaction <- bac_pup_ratio %>%
  #select(cohort, sex, dam_cage, diet, timepoint, treatment, sample) %>% # select all so can hash this out because all are selected by default.
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10")),
         cohort = factor(cohort),
         dam_cage = factor(dam_cage)
  ) %>%
  #group_by(metric) %>%
  nest() %>%
  mutate(lm = purrr::map(.x = data, ~ {
    if(var(.x$bac_bact_ratio) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(bac_bact_ratio ~ diet * treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(!data) %>%
  unnest(lm) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate))

write.csv(ratio_pups_diet_treatment_interaction, "ratio_pups_diet_treatment_interaction.csv", row.names = FALSE)

###Diet---
## Part 2: tests if diet (control or HFHS) has an effect in vehicle treatment group (no APC1472 or FOS+GOS) per ratio, sex, and timepoint.
diet_ratio_pups <- bac_pup_ratio %>%
  #select(cohort, dam_cage, sex, diet, timepoint, treatment, sample, Shannon, Simpson, Chao1) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(treatment == "Vehicle") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(diet_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$bac_bact_ratio) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(bac_bact_ratio ~ diet + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(diet_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "diet")

###Treatment
## Part 3: tests if treatment (vehicle, APC1472, or FOS+GOS) has an effect on HFHS diet group, per ratio, for sex and timepoint.
treatment_ratio_pups <- bac_pup_ratio %>%
  #select(cohort, dam_cage, sex, diet, timepoint, treatment, sample, Shannon, Simpson, Chao1) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(treatment = factor(x = treatment,
                            levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         diet = factor(x = diet,
                       levels = c("Control", "HFHS")),
         dam_cage = factor(dam_cage),
         cohort = factor(cohort), 
         sex = factor(sex, levels = c("Female", "Male")),
         timepoint = factor(timepoint, levels = c("week 5", "week 10"))
  ) %>%
  filter(diet == "HFHS") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(treatment_lmer = purrr::map(.x = data, .f = ~ {
    if(var(.x$bac_bact_ratio) == 0) {
      tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    } else {
      lmerTest::lmer(bac_bact_ratio ~ treatment + (1|cohort) + (1|dam_cage), data = .x) %>%
        broom.mixed::tidy()
    }
  })
  ) %>%
  select(-data) %>%  
  unnest(treatment_lmer) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  mutate(comparison_type = "treatment")

## Part 4: combining results from part 2 and 3
combined_diet_treatment_ratio_pups <- bind_rows(diet_ratio_pups, treatment_ratio_pups)

write.csv(combined_diet_treatment_ratio_pups,
          "combined_diet_treatment_ratio_pups.csv",
          row.names = FALSE)

###Restoration
## Part 5: Checks whether treatment actually reversed the effects of HFHS i.e restore microbiome back to control
diet_treatment_ratio_pups_restoration <- combined_diet_treatment_ratio_pups %>%
  filter(
    (comparison_type == "diet" & term == "dietHFHS") | 
      (comparison_type == "treatment" & term %in% c("treatmentFOS+GOS", "treatmentB. longum APC1472"))
  ) %>%
  group_by(sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[term == "dietHFHS"], default = NA),
    p.adj_diet = first(p.adj[term == "dietHFHS"], default = NA), 
    statistic_treatment_FOSGOS = first(statistic[term == "treatmentFOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[term == "treatmentFOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[term == "treatmentB. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[term == "treatmentB. longum APC1472"], default = NA)
  ) %>%
  mutate(
    restoration_FOSGOS = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_FOSGOS) &  
        p.adj_diet < 0.05 & p.adj_treatment_FOSGOS < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_FOSGOS),  
      "restored", "not restored"
    ),
    restoration_1472 = if_else(
      !is.na(p.adj_diet) & !is.na(p.adj_treatment_1472) &  
        p.adj_diet < 0.05 & p.adj_treatment_1472 < 0.05 &  
        sign(statistic_diet) != sign(statistic_treatment_1472),  
      "restored", "not restored"
    )
  )
write.csv(diet_treatment_ratio_pups_restoration, "diet_treatment_ratio_pups_restoration.csv", row.names = FALSE)
