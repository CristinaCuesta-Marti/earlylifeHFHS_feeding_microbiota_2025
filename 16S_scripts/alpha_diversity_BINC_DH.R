library(dada2)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(readxl)
library(tidyverse)
library(ggpubr)
library(vegan)
library(readxl)
library(writexl)

#metadata
metadata <- read_xlsx("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxee24/Phyloseq_maxee24_items/Metadata_BINC_maxee24_ps.xlsx")
metadata
meta <- column_to_rownames(metadata, var = "sample")
meta

#taxa table
tax <- read.csv("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxee24/Phyloseq_maxee24_items/taxa_maxee24_ps.csv", sep = ",", header = T, row.names = 1)
taxa <- as.matrix(tax)

#OTU table
seq <- read.csv("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxee24/Phyloseq_maxee24_items/seqtab_nochim_maxee24_ps.csv", header = TRUE, sep = ",", row.names = 1)
seqtab.no <- as.matrix(seq)

rownames(tax) == colnames(seq)

#make ps 
ps <- phyloseq(otu_table(seqtab.no, taxa_are_rows=FALSE),
               sample_data(meta),
               tax_table(taxa))
ps

##remove negblank from ps
pps <- subset_samples(ps, Sample_or_control != "control")
pps

#set diet colours
diet_colours <- c("Control Vehicle" = "#a0a0a4", "HFHS Vehicle" = "#f94040", "HFHS FOS+GOS" = "#addead", "HFHS B. longum APC1472" = "#e1c180")

#only bacteria
ps_bacteria <- pps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast" &
      Phylum != "Cyanobacteria/Chloroplast"
  )
ps_bacteria

##arrange diets, timepoints, sexes
#extract components from original phyloseq object
otu_table <- otu_table(ps_bacteria)
tax_table <- tax_table(ps_bacteria)
sample_d <- sample_data(ps_bacteria)
sample_d
sample_d <- data.frame(sample_d)
class(sample_d)
sample_d$timepoint <- factor(sample_d$timepoint, levels = c("w5", "w10"))
sample_d$sex <- factor(sample_d$sex, levels = c("M", "F"))
sample_d$diet <- factor(sample_d$diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS", "HFHS B. longum APC1472"))
str(sample_d)
#convert the modified sample data back to sample_data
sample_d <- sample_data(sample_d)
# Create a new phyloseq object with the reordered sample data
ps_reordered <- phyloseq(otu_table, tax_table, sample_d)

################## Dams
# filter for dams only
ps_dam <- subset_samples(ps_reordered, age != "pup")

## Simpson
simp_dam <- plot_richness(ps_dam, x="diet", measures=c("Simpson")) + 
  geom_boxplot(aes(fill = diet)) +
  scale_fill_manual(values = diet_colours) +
  theme(strip.background = element_rect(fill = "white")) +
  geom_pwc(label = "p.format", p.adjust.method = "BH",
           #label = "p.adj.signif",
           step.increase = 0.07, vjust = 0.3, hide.ns = T) +
  theme_bw() + ggtitle("Simpson diversity dams ") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size=10),
        axis.text = element_text(size = 10), legend.text = element_text(size=10),
        axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

## Shannon
shan_dam <- plot_richness(ps_dam, x="diet", measures=c("Shannon")) + 
  geom_boxplot(aes(fill = diet)) +
  scale_fill_manual(values = diet_colours) +
  theme(strip.background = element_rect(fill = "white")) +
  geom_pwc(label = "p.format", p.adjust.method = "BH", 
           step.increase = 0.07, vjust = 0.3, hide.ns = T) +
  theme_bw() + ggtitle("Shannon diversity dams ") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size=10),
        axis.text = element_text(size = 10), legend.text = element_text(size=10),
        axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

#?geom_pwc()

## Chao1
chao_dam <- plot_richness(ps_dam, x="diet", measures=c("Chao1")) + 
  geom_boxplot(aes(fill = diet)) +
  scale_fill_manual(values = diet_colours) +
  theme(strip.background = element_rect(fill = "white")) +
  geom_pwc(label = "p.format", p.adjust.method = "BH", 
           step.increase = 0.07, vjust = 0.3, hide.ns = T) +
  theme_bw() + ggtitle("Chao1 diversity dams ") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size=10),
        axis.text = element_text(size = 10), legend.text = element_text(size=10),
        axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

#setwd("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Alpha Diversity")
ggsave("Simpson_dams.svg", simp_dam, width = 12, height = 8, device = "svg")
ggsave("Simpson_dams.png", simp_dam, width = 12, height = 8, dpi = 1200)

ggsave("Shannon_dams.svg", shan_dam, width = 12, height = 8, device = "svg")
ggsave("Shannon_dams.png", shan_dam, width = 12, height = 8, dpi = 1200)

ggsave("Chao1_dams.svg", chao_dam, width = 12, height = 8, device = "svg")
ggsave("Chao1_dams.png", chao_dam, width = 12, height = 8, dpi = 1200)

################## Pups
## I will only make plots timepoint separated. Male and female beside each other. Diet on x axis

######## Pups w5
# filter for pups w5 only
ps_pupw5 <- subset_samples(ps_reordered, timepoint == "w5")

## Simpson
simp_pupw5 <- plot_richness(ps_pupw5, x="diet", measures=c("Simpson")) + 
  geom_boxplot(aes(fill = diet)) +
  facet_grid(~ sex) +
  scale_fill_manual(values = diet_colours) +
  theme(strip.background = element_rect(fill = "white")) +
  geom_pwc(label = "p.format", p.adjust.method = "BH",
           #label = "p.adj.signif",
           step.increase = 0.07, vjust = 0.3, hide.ns = T) +
  theme_bw() + ggtitle("Simpson diversity pups week 5 ") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size=10),
        axis.text = element_text(size = 10), legend.text = element_text(size=10),
        axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

## Shannon
shan_pupw5 <- plot_richness(ps_pupw5, x="diet", measures=c("Shannon")) + 
  geom_boxplot(aes(fill = diet)) +
  facet_grid(~ sex) +
  scale_fill_manual(values = diet_colours) +
  theme(strip.background = element_rect(fill = "white")) +
  geom_pwc(label = "p.format", p.adjust.method = "BH",
           #label = "p.adj.signif",
           step.increase = 0.07, vjust = 0.3, hide.ns = T) +
  theme_bw() + ggtitle("Shannon diversity pups week 5 ") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size=10),
        axis.text = element_text(size = 10), legend.text = element_text(size=10),
        axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

## Chao1
chao_pupw5 <- plot_richness(ps_pupw5, x="diet", measures=c("Chao1")) + 
  geom_boxplot(aes(fill = diet)) +
  facet_grid(~ sex) +
  scale_fill_manual(values = diet_colours) +
  theme(strip.background = element_rect(fill = "white")) +
  geom_pwc(label = "p.format", p.adjust.method = "BH",
           #label = "p.adj.signif",
           step.increase = 0.07, vjust = 0.3, hide.ns = T) +
  theme_bw() + ggtitle("Chao1 diversity pups week 5 ") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size=10),
        axis.text = element_text(size = 10), legend.text = element_text(size=10),
        axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

#setwd("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Alpha Diversity")
ggsave("Simpson_pups_w5.svg", simp_pupw5, width = 12, height = 8, device = "svg")
ggsave("Simpson_pups_w5.png", simp_pupw5, width = 12, height = 8, dpi = 1200)

ggsave("Shannon_pups_w5.svg", shan_pupw5, width = 12, height = 8, device = "svg")
ggsave("Shannon_pups_w5.png", shan_pupw5, width = 12, height = 8, dpi = 1200)

ggsave("Chao1_pups_w5.svg", chao_pupw5, width = 12, height = 8, device = "svg")
ggsave("Chao1_pups_w5.png", chao_pupw5, width = 12, height = 8, dpi = 1200)

######## Pups w10
# filter for pups w10 only
ps_pupw10 <- subset_samples(ps_reordered, timepoint == "w10")

## Simpson
simp_pupw10 <- plot_richness(ps_pupw10, x="diet", measures=c("Simpson")) + 
  geom_boxplot(aes(fill = diet)) +
  facet_grid(~ sex) +
  scale_fill_manual(values = diet_colours) +
  theme(strip.background = element_rect(fill = "white")) +
  geom_pwc(label = "p.format", p.adjust.method = "BH",
           #label = "p.adj.signif",
           step.increase = 0.07, vjust = 0.3, hide.ns = T) +
  theme_bw() + ggtitle("Simpson diversity pups week 10 ") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size=10),
        axis.text = element_text(size = 10), legend.text = element_text(size=10),
        axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

## Shannon
shan_pupw10 <- plot_richness(ps_pupw10, x="diet", measures=c("Shannon")) + 
  geom_boxplot(aes(fill = diet)) +
  facet_grid(~ sex) +
  scale_fill_manual(values = diet_colours) +
  theme(strip.background = element_rect(fill = "white")) +
  geom_pwc(label = "p.format", p.adjust.method = "BH",
           #label = "p.adj.signif",
           step.increase = 0.07, vjust = 0.3, hide.ns = T) +
  theme_bw() + ggtitle("Shannon diversity pups week 10 ") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size=10),
        axis.text = element_text(size = 10), legend.text = element_text(size=10),
        axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

## Chao1
chao_pupw10 <- plot_richness(ps_pupw10, x="diet", measures=c("Chao1")) + 
  geom_boxplot(aes(fill = diet)) +
  facet_grid(~ sex) +
  scale_fill_manual(values = diet_colours) +
  theme(strip.background = element_rect(fill = "white")) +
  geom_pwc(label = "p.format", p.adjust.method = "BH",
           #label = "p.adj.signif",
           step.increase = 0.07, vjust = 0.3, hide.ns = T) +
  theme_bw() + ggtitle("Chao1 diversity pups week 10 ") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size=10),
        axis.text = element_text(size = 10), legend.text = element_text(size=10),
        axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

#setwd("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Alpha Diversity")
ggsave("Simpson_pups_w10.svg", simp_pupw10, width = 12, height = 8, device = "svg")
ggsave("Simpson_pups_w10.png", simp_pupw10, width = 12, height = 8, dpi = 1200)

ggsave("Shannon_pups_w10.svg", shan_pupw10, width = 12, height = 8, device = "svg")
ggsave("Shannon_pups_w10.png", shan_pupw10, width = 12, height = 8, dpi = 1200)

ggsave("Chao1_pups_w10.svg", chao_pupw10, width = 12, height = 8, device = "svg")
ggsave("Chao1_pups_w10.png", chao_pupw10, width = 12, height = 8, dpi = 1200)