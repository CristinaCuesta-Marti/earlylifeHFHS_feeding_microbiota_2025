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
#R.version.string
#BiocManager::install("phyloseq")

###Making the phyloseq object for downstream analysis###
setwd("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items")
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

#-----------------------------------------------------------------------------------------
#####Beta diversity pups
diet_colours <- c("Control Vehicle" = "#a0a0a4", "HFHS Vehicle" = "#f94040", "HFHS FOS+GOS" = "#addead", "HFHS B. longum APC1472" = "#e1c180")
#Extract components from the original phyloseq object
otu_table <- otu_table(PSr_pup)
tax_table <- tax_table(PSr_pup)
# Extract and modify the sample data to include the reordered diet levels
sample_d <- sample_data(PSr_pup)
sample_d
sample_d <- as.data.frame(sample_d)
sample_d$diet <- factor(sample_d$diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS", "HFHS B. longum APC1472"))
sample_d$timepoint <- factor(sample_d$timepoint, levels = c("w5", "w10"))
sample_d$sex <- factor(sample_d$sex, levels = c("M", "F"))
PSr_pup_reordered <- phyloseq(otu_table, tax_table, sample_d)
sample_data(PSr_pup_reordered)
class(PSr_pup_reordered)


####subsetting w5 male
PSr_pup_reordered_w5M <- subset_samples(PSr_pup_reordered, sample_data(PSr_pup_reordered)$timepoint !="w10" 
                                       & sample_data(PSr_pup_reordered)$sex != "F")
####subsetting w5 female
PSr_pup_reordered_w5F <- subset_samples(PSr_pup_reordered, sample_data(PSr_pup_reordered)$timepoint !="w10" 
                                        & sample_data(PSr_pup_reordered)$sex != "M")
####subsetting w10 male
PSr_pup_reordered_w10M <- subset_samples(PSr_pup_reordered, sample_data(PSr_pup_reordered)$timepoint !="w5" 
                                        & sample_data(PSr_pup_reordered)$sex != "F")
####subsetting w10 female
PSr_pup_reordered_w10F <- subset_samples(PSr_pup_reordered, sample_data(PSr_pup_reordered)$timepoint !="w5" 
                                        & sample_data(PSr_pup_reordered)$sex != "M")


ordu.bray_w5M = ordinate(PSr_pup_reordered_w5M, "PCoA", "bray")
scree.plot_w5M <- plot_scree(ordu.bray_w5M, "Check for importance of axis in Scree plot W5 MALE")
print(scree.plot_w5M)

ordu.bray_w5F = ordinate(PSr_pup_reordered_w5F, "PCoA", "bray")
scree.plot_w5F <- plot_scree(ordu.bray_w5F, "Check for importance of axis in Scree plot W5 FEMALE")
print(scree.plot_w5F)

ordu.bray_w10M = ordinate(PSr_pup_reordered_w10M, "PCoA", "bray")
scree.plot_w10M <- plot_scree(ordu.bray_w10M, "Check for importance of axis in Scree plot W10 MALE")
print(scree.plot_w10M)

ordu.bray_w10F = ordinate(PSr_pup_reordered_w10F, "PCoA", "bray")
scree.plot_w10F <- plot_scree(ordu.bray_w10F, "Check for importance of axis in Scree plot W10 FEMALE")
print(scree.plot_w10F)

setwd("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Beta_diversity_calculatedrelA/Beta Diversity 2025")
#install.packages("svglite")
getwd()
ggsave("scree_pups_w5M.svg", scree.plot_w5M, width = 12, height = 8, device = "svg")
ggsave("scree_pups_w5M.png", scree.plot_w5M, width = 12, height = 8, dpi = 1200)
ggsave("scree_pups_w5F.svg", scree.plot_w5F, width = 12, height = 8, device = "svg")
ggsave("scree_pups_w5F.png", scree.plot_w5F, width = 12, height = 8, dpi = 1200)
ggsave("scree_pups_w10M.svg", scree.plot_w10M, width = 12, height = 8, device = "svg")
ggsave("scree_pups_w10M.png", scree.plot_w10M, width = 12, height = 8, dpi = 1200)
ggsave("scree_pups_w10F.svg", scree.plot_w10F, width = 12, height = 8, device = "svg")
ggsave("scree_pups_w10F.png", scree.plot_w10F, width = 12, height = 8, dpi = 1200)

#PCoA for w5 male 
abund_w5M <- otu_table(PSr_pup_reordered_w5M)
abund_ranks_w5M <- t(apply(abund_w5M, 1, rank))
ps_ranks1_w5M <- plot_ordination(PSr_pup_reordered_w5M, ordinate(PSr_pup_reordered_w5M, "PCoA", "bray"), color = "diet") + 
  geom_point(size = 1) + 
  stat_ellipse(aes(group = diet), linewidth = 1) +
  scale_color_manual(values = diet_colours)
ps_ranks1_w5M

bray_pcoa_w5M <-  ps_ranks1_w5M + 
  theme_bw() +
  ggtitle("PCoA plot Bray-Curtis - pups w5 Male") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11), 
        axis.text = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        panel.border = element_blank(),
        panel.background = element_blank(),  # Remove the panel background shading
        strip.background = element_blank())  # Remove facet panel background shading
  #facet_grid(sex ~ timepoint)
bray_pcoa_w5M
ggsave("beta_bray_w5M.svg", bray_pcoa_w5M, width = 12, height = 8, device = "svg")
ggsave("beta_bray_w5M.png", bray_pcoa_w5M, width = 12, height = 8, dpi = 1200)

#PCoA for w5 female
abund_w5F <- otu_table(PSr_pup_reordered_w5F)
abund_ranks_w5F <- t(apply(abund_w5F, 1, rank))
ps_ranks1_w5F <- plot_ordination(PSr_pup_reordered_w5F, ordinate(PSr_pup_reordered_w5F, "PCoA", "bray"), color = "diet") + 
  geom_point(size = 1) + 
  stat_ellipse(aes(group = diet), linewidth = 1) +
  scale_color_manual(values = diet_colours)
ps_ranks1_w5F

bray_pcoa_w5F <-  ps_ranks1_w5F + 
  theme_bw() +
  ggtitle("PCoA plot Bray-Curtis - pups w5 Female") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11), 
        axis.text = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        panel.border = element_blank(),
        panel.background = element_blank(),  # Remove the panel background shading
        strip.background = element_blank())  # Remove facet panel background shading
#facet_grid(sex ~ timepoint)
bray_pcoa_w5F
ggsave("beta_bray_w5F.svg", bray_pcoa_w5F, width = 12, height = 8, device = "svg")
ggsave("beta_bray_w5F.png", bray_pcoa_w5F, width = 12, height = 8, dpi = 1200)

#PCoA for w10 male 
abund_w10M <- otu_table(PSr_pup_reordered_w10M)
abund_ranks_w10M <- t(apply(abund_w10M, 1, rank))
ps_ranks1_w10M <- plot_ordination(PSr_pup_reordered_w10M, ordinate(PSr_pup_reordered_w10M, "PCoA", "bray"), color = "diet") + 
  geom_point(size = 1) + 
  stat_ellipse(aes(group = diet), linewidth = 1) +
  scale_color_manual(values = diet_colours)
ps_ranks1_w10M

bray_pcoa_w10M <-  ps_ranks1_w10M + 
  theme_bw() +
  ggtitle("PCoA plot Bray-Curtis - pups w10 Male") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11), 
        axis.text = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        panel.border = element_blank(),
        panel.background = element_blank(),  # Remove the panel background shading
        strip.background = element_blank())  # Remove facet panel background shading
#facet_grid(sex ~ timepoint)
bray_pcoa_w10M
ggsave("beta_bray_w10M.svg", bray_pcoa_w10M, width = 12, height = 8, device = "svg")
ggsave("beta_bray_w10M.png", bray_pcoa_w10M, width = 12, height = 8, dpi = 1200)

#PCoA for w10 female 
abund_w10F <- otu_table(PSr_pup_reordered_w10F)
abund_ranks_w10F <- t(apply(abund_w10F, 1, rank))
ps_ranks1_w10F <- plot_ordination(PSr_pup_reordered_w10F, ordinate(PSr_pup_reordered_w10F, "PCoA", "bray"), color = "diet") + 
  geom_point(size = 1) + 
  stat_ellipse(aes(group = diet), linewidth = 1) +
  scale_color_manual(values = diet_colours)
ps_ranks1_w10F

bray_pcoa_w10F <-  ps_ranks1_w10F + 
  theme_bw() +
  ggtitle("PCoA plot Bray-Curtis - pups w10 Female") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11), 
        axis.text = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        panel.border = element_blank(),
        panel.background = element_blank(),  # Remove the panel background shading
        strip.background = element_blank())  # Remove facet panel background shading
#facet_grid(sex ~ timepoint)
bray_pcoa_w10F
ggsave("beta_bray_w10F.svg", bray_pcoa_w10F, width = 12, height = 8, device = "svg")
ggsave("beta_bray_w10F.png", bray_pcoa_w10F, width = 12, height = 8, dpi = 1200)

#####Put the four PCoAs on the same plot with separate scales
#or not










#########################Dams 
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

ordu.bray2d = ordinate(PSr_dam_reordered, "PCoA", "bray")
scree.plotd <- plot_scree(ordu.bray2d, "Check for importance of axis in Scree plot DAMS")
print(scree.plotd)
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Calculated_relative_abundances")
#install.packages("svglite")
ggsave("scree_dams_new.svg", scree.plotd, width = 12, height = 8, device = "svg")
ggsave("scree_dams_new.png", scree.plotd, width = 12, height = 8, dpi = 1200)

abund_alld <- otu_table(PSr_dam_reordered)
abund_ranks_alld <- t(apply(abund_alld, 1, rank))
ps_ranks1_alld <- plot_ordination(PSr_dam_reordered, ordinate(PSr_dam_reordered, "PCoA", "bray"), color = "diet") + 
  geom_point(size = 3) + 
  stat_ellipse(aes(group = diet), linewidth = 1) +
  scale_color_manual(values = diet_colours)
ps_ranks1_alld
diet_colours <- c("Control Vehicle" = "#a0a0a4", "HFHS Vehicle" = "#f94040", "HFHS FOS+GOS" = "#addead", "HFHS B. longum APC1472" = "#e1c180")

bray_pcoa_dams <-  ps_ranks1_alld + 
  theme_bw() +
  ggtitle("PCoA plot Bray-Curtis - dams") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11), 
        axis.text = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        panel.border = element_blank(),
        panel.background = element_blank(),  # Remove the panel background shading
        strip.background = element_blank())   # Remove facet panel background shading
#facet_grid(sex ~ timepoint)
bray_pcoa_dams
ggsave("beta_bray_dams.svg", bray_pcoa_dams, width = 12, height = 8, device = "svg")
ggsave("beta_bray_dams.png", bray_pcoa_dams, width = 12, height = 8, dpi = 1200)

###Script from Dhrati. 
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", force = T) 
library(devtools)
library(remotes)
library(pairwiseAdonis)


