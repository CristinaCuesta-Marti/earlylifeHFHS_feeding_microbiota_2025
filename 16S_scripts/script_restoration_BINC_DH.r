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
#install.packages("janitor")
library(janitor)
library(tibble)
#install.packages("purrr")
library(purrr)
library(rstatix)
library(data.table)

###make a phyloseq object 

ps3

##separate the pups and dams, and make do separate ps based on age column

ps_pup <- subset_samples(ps3, sample_data(ps)$age =="pup")
ps_dam <- subset_samples(ps3, sample_data(ps)$age =="dam")


#Make relative abundances
PSr_pup <- transform_sample_counts(ps_pup, function(x) x/sum(x))
PSr_dam <- transform_sample_counts(ps_dam, function(x) x/sum(x))

########################### Alpha Diversity
##Calculate alpha diversity
alpha_diversity <- estimate_richness(ps_pup, measures = c("Shannon", "Chao1", "Simpson")) # estimate_richness needs a phyloseq
rownames(alpha_diversity) <- gsub("^X", "", rownames(alpha_diversity)) # remove pesky X

# Merge with sample metadata VERY IMPORTANT STEP
meta <- sample_data(ps_pup)
alpha_diversity <- merge(meta, alpha_diversity, by = "row.names")
row.names(meta)
row.names(alpha_diversity)

#lm_model <- lm(Chao1 ~ diet * treatment * sex * timepoint, data = alpha_diversity)
#summary(lm_model)
### The above script may work. I'm going to stick it into Cristina's script now. I get the same outputs Cristina and my script

########## Interactions section 1
#Shannon
alpha_diversity_interactions <- alpha_diversity %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Shannon) %>%
  as_tibble() %>%
  clean_names() %>%
  #group_by() %>%
  nest() %>%
  mutate(lm = map(.x = data,
                         .f = ~ lm(shannon ~ diet * treatment * sex * timepoint,
                                   data = .x) %>%
                           broom::tidy())) %>%
  select(!data) %>%
  unnest(lm) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  write.csv(., "Alpha_diversity_shannon_interactions.csv", row.names = FALSE)

#setwd("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Cristina's restoration script 2025 - outputs only")

#Simpson
alpha_diversity_interactions <- alpha_diversity %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Simpson) %>%
  as_tibble() %>%
  clean_names() %>%
  #group_by() %>%
  nest() %>%
  mutate(lm = map(.x = data,
                  .f = ~ lm(simpson ~ diet * treatment * sex * timepoint,
                            data = .x) %>%
                    broom::tidy())) %>%
  select(!data) %>%
  unnest(lm) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  write.csv(., "Alpha_diversity_simpson_interactions.csv", row.names = FALSE)

#Chao1
alpha_diversity_interactions <- alpha_diversity %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Chao1) %>%
  as_tibble() %>%
  clean_names() %>%
  #group_by() %>%
  nest() %>%
  mutate(lm = map(.x = data,
                  .f = ~ lm(chao1 ~ diet * treatment * sex * timepoint,
                            data = .x) %>%
                    broom::tidy())) %>%
  select(!data) %>%
  unnest(lm) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) %>%
  write.csv(., "Alpha_diversity_chao1_interactions.csv", row.names = FALSE)

####### Pairwise differences between diet groups and treatments. Analysed separately
#Shannon
diet_alpha_diversity_shan <- alpha_diversity %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Shannon) %>%
  as_tibble() %>%
  clean_names() %>%
  filter(treatment == "Vehicle") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(diet_alpha = map(.x = data,
                          .f = ~ rstatix::pairwise_t_test(data = .x,
                                                          formula = shannon ~ diet,
                                                          pool.sd = FALSE))) %>%
  unnest(c(diet_alpha)) %>%
  select(-data) %>%
  mutate(comparison_type = "diet")

treatment_alpha_diversity_shan <- alpha_diversity %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Shannon) %>%
  as_tibble() %>%
  clean_names() %>%
  filter(diet == "HFHS") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(treatment_alpha = map(.x = data,
                          .f = ~ rstatix::pairwise_t_test(data = .x,
                                                          formula = shannon ~ treatment,
                                                          pool.sd = FALSE))) %>%
  unnest(c(treatment_alpha)) %>%
  select(-data) %>%
  mutate(comparison_type = "treatment")

combined_diet_treatment_alpha_shan <- rbind(diet_alpha_diversity_shan, treatment_alpha_diversity_shan)
write.csv(combined_diet_treatment_alpha_shan, "combined_diet_treatment_alpha_shannon.csv", row.names = FALSE)

#Simpson
diet_alpha_diversity_simp <- alpha_diversity %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Simpson) %>%
  as_tibble() %>%
  clean_names() %>%
  filter(treatment == "Vehicle") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(diet_alpha = map(.x = data,
                          .f = ~ rstatix::pairwise_t_test(data = .x,
                                                          formula = simpson ~ diet,
                                                          pool.sd = FALSE))) %>%
  unnest(c(diet_alpha)) %>%
  select(-data) %>%
  mutate(comparison_type = "diet")

treatment_alpha_diversity_simp <- alpha_diversity %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Simpson) %>%
  as_tibble() %>%
  clean_names() %>%
  filter(diet == "HFHS") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(treatment_alpha = map(.x = data,
                               .f = ~ rstatix::pairwise_t_test(data = .x,
                                                               formula = simpson ~ treatment,
                                                               pool.sd = FALSE))) %>%
  unnest(c(treatment_alpha)) %>%
  select(-data) %>%
  mutate(comparison_type = "treatment")

combined_diet_treatment_alpha_simp <- rbind(diet_alpha_diversity_simp, treatment_alpha_diversity_simp)
write.csv(combined_diet_treatment_alpha_simp, "combined_diet_treatment_alpha_simpson.csv", row.names = FALSE)

#Chao1
diet_alpha_diversity_chao <- alpha_diversity %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Chao1) %>%
  as_tibble() %>%
  clean_names() %>%
  filter(treatment == "Vehicle") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(diet_alpha = map(.x = data,
                          .f = ~ rstatix::pairwise_t_test(data = .x,
                                                          formula = chao1 ~ diet,
                                                          pool.sd = FALSE))) %>%
  unnest(c(diet_alpha)) %>%
  select(-data) %>%
  mutate(comparison_type = "diet")

treatment_alpha_diversity_chao <- alpha_diversity %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Chao1) %>%
  as_tibble() %>%
  clean_names() %>%
  filter(diet == "HFHS") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(treatment_alpha = map(.x = data,
                               .f = ~ rstatix::pairwise_t_test(data = .x,
                                                               formula = chao1 ~ treatment,
                                                               pool.sd = FALSE))) %>%
  unnest(c(treatment_alpha)) %>%
  select(-data) %>%
  mutate(comparison_type = "treatment")

combined_diet_treatment_alpha_chao <- rbind(diet_alpha_diversity_chao, treatment_alpha_diversity_chao)
write.csv(combined_diet_treatment_alpha_chao, "combined_diet_treatment_alpha_chao1.csv", row.names = FALSE)

######## Section 3 which checks if either treatment restores a difference brought about by HFHS.
#Shannon
restore_diet_treatment_alpha_shan <- combined_diet_treatment_alpha_shan %>%
  filter(
    (comparison_type == "diet" & group1 == "Control" & group2 == "HFHS") |
  (comparison_type == "treatment" & group1 == "Vehicle" & group2 %in% c("FOS+GOS", "B. longum APC1472"))
  ) %>%
  group_by(sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[comparison_type == "diet"], default = NA),
    p.adj_diet = first(p.adj[comparison_type == "diet"], default = NA),
    statistic_treatment_FOSGOS = first(statistic[comparison_type == "treatment" & group2 == "FOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[comparison_type == "treatment" & group2 == "FOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[comparison_type == "treatment" & group2 == "B. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[comparison_type == "treatment" & group2 == "B. longum APC1472"], default = NA)
  ) %>%
  mutate(restoration_FOSGOS = if_else(
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

write.csv(restore_diet_treatment_alpha_shan, "restore_diet_treatment_alpha_shannon.csv", row.names = FALSE)

#Simpson
restore_diet_treatment_alpha_simp <- combined_diet_treatment_alpha_simp %>%
  filter(
    (comparison_type == "diet" & group1 == "Control" & group2 == "HFHS") |
      (comparison_type == "treatment" & group1 == "Vehicle" & group2 %in% c("FOS+GOS", "B. longum APC1472"))
  ) %>%
  group_by(sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[comparison_type == "diet"], default = NA),
    p.adj_diet = first(p.adj[comparison_type == "diet"], default = NA),
    statistic_treatment_FOSGOS = first(statistic[comparison_type == "treatment" & group2 == "FOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[comparison_type == "treatment" & group2 == "FOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[comparison_type == "treatment" & group2 == "B. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[comparison_type == "treatment" & group2 == "B. longum APC1472"], default = NA)
  ) %>%
  mutate(restoration_FOSGOS = if_else(
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

write.csv(restore_diet_treatment_alpha_simp, "restore_diet_treatment_alpha_simpson.csv", row.names = FALSE)

#Chao1
restore_diet_treatment_alpha_chao <- combined_diet_treatment_alpha_chao %>%
  filter(
    (comparison_type == "diet" & group1 == "Control" & group2 == "HFHS") |
      (comparison_type == "treatment" & group1 == "Vehicle" & group2 %in% c("FOS+GOS", "B. longum APC1472"))
  ) %>%
  group_by(sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[comparison_type == "diet"], default = NA),
    p.adj_diet = first(p.adj[comparison_type == "diet"], default = NA),
    statistic_treatment_FOSGOS = first(statistic[comparison_type == "treatment" & group2 == "FOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[comparison_type == "treatment" & group2 == "FOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[comparison_type == "treatment" & group2 == "B. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[comparison_type == "treatment" & group2 == "B. longum APC1472"], default = NA)
  ) %>%
  mutate(restoration_FOSGOS = if_else(
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

write.csv(restore_diet_treatment_alpha_chao, "restore_diet_treatment_alpha_chao1.csv", row.names = FALSE)

##################################### Beta diversity
##Calculate beta diversity distance

#### Merge with sample metadata VERY IMPORTANT STEP
#Calculate Bray
DistBC <- distance(PSr_pup, method = "bray")

#Perform PCoA on each distance matrix
ordBC <- ordinate(PSr_pup, method = "PCoA", distance = DistBC)

#Check variance of each axis
ordBC$values$Eigenvalues
eig_values <- ordBC$values$Eigenvalues
percent_explained <- (eig_values / sum(eig_values)) * 100
percent_explained[1:5] # Axis 1 accounts for 25% of variance
#plot_scree(ordBC, "Scree Plot: Bray-Curtis PSr_pup")

#Let's do some plotting
#plot_ordination(PSr_pup, ordBC, color = "diet") + 
 # geom_point() +
  #ggtitle("PCoA: Bray-Curtis")

#convert phyloseq sample data dataframe to a data table
sdt <- data.table(as(sample_data(PSr_pup), "data.frame"),
                  keep.rownames = TRUE) # rn (rownames) matches the 196 sample names
setnames(sdt, "rn", "SampleID") # rn is now called SampleID
                  
##Join sample data and ordination axes together in one data.table
ordBCdt <- data.table(as.data.frame(ordBC$vectors), keep.rownames = TRUE)
setnames(ordBCdt, "rn", "SampleID")
#str(ordBCdt)
#head(ordBCdt)
#colnames(ordBCdt)
#colnames(sdt)
#all(ordBCdt$SampleID %in% sdt$SampleID)
#all(sdt$SampleID %in% ordBCdt$SampleID)
setkey(ordBCdt, SampleID)
setkey(sdt, SampleID)
#ordBCsdt <- ordBCdt[sdt]
ordBCsdt <- sdt[ordBCdt]
#setorder(ordBCsdt, timepoint)
str(ordBCsdt)

### Section 1 interactions script beta diversity. Only do for axes 1.
## Axis 1 captures the most variation (25%) in my dataset, followed by axis 2 and so on.
beta_diversity_analysis <- ordBCsdt %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Axis.1, Axis.2, Axis.3) %>%
  as_tibble() %>%
  clean_names() %>%
  nest() %>%
  mutate(lm = map(.x = data,
                  .f = ~ lm(axis_1 # use axis_1 instead of Axis.1 because clean_names
                            ~ diet * treatment * sex * timepoint, 
                            data = .x) %>%
                    broom::tidy())) %>%
  select(!data) %>%
  unnest(lm) %>%
  filter(term != "(Intercept)" & !is.na(estimate)) 
  #write.csv(., "beta_diversity_interactions.csv", row.names = FALSE)

getwd()
#setwd("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Cristina's restoration script 2025 - outputs only")

### Section 2 pairwise differences between diets and treatments
diet_beta_diversity <- ordBCsdt %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Axis.1) %>%
  as_tibble() %>%
  clean_names() %>%
  filter(treatment == "Vehicle") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(diet_beta = map(.x = data,
                         .f = ~ pairwise_wilcox_test(data = .x,
                                                                   formula = axis_1 ~ diet,
                                                                   p.adjust.method = "BH"))) %>%
  unnest(c(diet_beta)) %>%
  select(-data) %>%
  mutate(comparison_type = "diet")

treatment_beta_diversity <- ordBCsdt %>%
  select(sex, cohort, diet, treatment, timepoint, sample, Axis.1) %>%
  as_tibble() %>%
  clean_names() %>%
  filter(diet == "HFHS") %>%
  group_by(sex, timepoint) %>%
  nest() %>%
  mutate(treatment_beta = map(.x = data,
                              .f = ~ pairwise_wilcox_test(data = .x,
                                                                        formula = axis_1 ~ treatment,
                                                                        p.adjust.method = "BH"))) %>%
  unnest(c(treatment_beta)) %>%
  select(-data) %>%
  mutate(comparison_type = "treatment")

combined_diet_treatment_beta <- rbind(diet_beta_diversity, treatment_beta_diversity)
write.csv(combined_diet_treatment_beta, "combined_diet_treatment_beta.csv", row.names = FALSE)

###Section 3 script for restoration effects
restore_diet_treatment_beta <- combined_diet_treatment_beta %>%
  filter(
    (comparison_type == "diet" & group1 == "Control" & group2 == "HFHS") |
      (comparison_type == "treatment" & group1 == "Vehicle" & group2 %in% c("FOS+GOS", "B. longum APC1472"))
  ) %>%
  group_by(sex, timepoint) %>%
  summarise(
    statistic_diet = first(statistic[comparison_type == "diet"], default = NA),
    p.adj_diet = first(p.adj[comparison_type == "diet"], default = NA),
    statistic_treatment_FOSGOS = first(statistic[comparison_type == "treatment" & group2 == "FOS+GOS"], default = NA),
    p.adj_treatment_FOSGOS = first(p.adj[comparison_type == "treatment" & group2 == "FOS+GOS"], default = NA),
    statistic_treatment_1472 = first(statistic[comparison_type == "treatment" & group2 == "B. longum APC1472"], default = NA),
    p.adj_treatment_1472 = first(p.adj[comparison_type == "treatment" & group2 == "B. longum APC1472"], default = NA)
  ) %>%
  mutate(restoration_FOSGOS = if_else(
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
write.csv(restore_diet_treatment_beta, "restore_diet_treatment_beta.csv", row.names = FALSE)

