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
library(writexl)
library(svglite)
library(rstatix)

setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items")

#Metadata or sample_data
samdf <- read_excel("Metadata_BINC_maxee24_ps.xlsx")
sam2 <- samdf %>% remove_rownames %>% column_to_rownames(var="sample")

#taxa table or tax_table
taxa <- read.csv("taxa_maxee24_ps.csv", sep = ",", row.names = 1)
taxa_mat <- as.matrix(taxa)

#sequence table or otu_table
seqtab.nochim <- read.csv("seqtab_nochim_maxee24_ps.csv", sep = ",", row.names = 1)
seqtab.nochim_t <- t(seqtab.nochim)
head(colnames(seqtab.nochim))
head(rownames(seqtab.nochim))
seq_mat <- as.matrix(seqtab.nochim)

#make phyloseq object
ps <- phyloseq(otu_table(seq_mat, taxa_are_rows=FALSE), 
               sample_data(sam2), 
               tax_table(taxa_mat))
ps

##I think Dara had nothing in decontam, so removing blanks
##thus removing those for now
ps.noblank <- subset_samples(ps, sample_data(ps)$Sample_or_control !="control")

#remove non-bacteria

ps3 <- ps.noblank %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast" &
      Phylum != "Cyanobacteria/Chloroplast"
  )

ps3

#arrange diets
str(ps3)
#Extract components from the original phyloseq object
otu_table <- otu_table(ps3)
tax_table <- tax_table(ps3)
# Extract and modify the sample data to include the reordered diet levels
sample_data_reordered <- sample_data(ps3)
sample_data_reordered <- as.data.frame(sample_data_reordered)
sample_data_reordered$diet <- factor(sample_data_reordered$diet, levels = c("Control Vehicle", "HFHS Vehicle"
                                                                            , "HFHS FOS+GOS", "HFHS B. longum APC1472"))
sample_data_reordered$timepoint <- factor(sample_data_reordered$timepoint, levels = c("w5", "w10"))
sample_data_reordered$sex <- factor(sample_data_reordered$sex, levels = c("M", "F"))
sample_data <- sample_data_reordered
ps3_reordered <- phyloseq(otu_table, tax_table, sample_data)
sample_data(ps3_reordered)
class(ps3_reordered)

##separate the pups and dams, and make do separate ps based on age column

ps_pup <- subset_samples(ps3_reordered, sample_data(ps3_reordered)$age =="pup")
ps_dam <- subset_samples(ps3_reordered, sample_data(ps3_reordered)$age =="dam")

#Make relative abundances
PSr_pup <- transform_sample_counts(ps_pup, function(x) x/sum(x))
PSr_dam <- transform_sample_counts(ps_dam, function(x) x/sum(x))

sample_sums <- sample_sums(PSr_pup)
all(sample_sums == 1)

#########################pups
#####Phylum
###Bars
ps_pup_phylum <- PSr_pup %>%
  tax_glom(taxrank = "Phylum") %>%                    
  psmelt() %>%                                         
  #filter(Abundance > 0.02) %>%               
  arrange(Phylum) 
str(ps_pup_phylum)

ps_pup_phylum$Phylum <- factor(ps_pup_phylum$Phylum, levels = ps_pup_phylum %>%
                                 group_by(Phylum) %>%
                                 summarise(total_abundance = sum(Abundance)) %>%
                                 arrange(desc(total_abundance)) %>%
                                 pull(Phylum))

###heatmap checks
ps_pup_phylum_z <- ps_pup_phylum %>%
  group_by(Phylum, sex, timepoint, diet) %>%
  mutate(mean = mean(Abundance)) %>%
  mutate(SD = sd(Abundance)) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

sum(is.na(ps_pup_phylum_z$z_score))

OTU_actinobacteriota <- ps_pup_phylum_z %>%
  filter(Phylum == "Actinobacteriota", sex == "M", timepoint == "w5", diet == "Control Vehicle")

OTU_deferribacterota <- ps_pup_phylum_z %>%
  filter(Phylum == "Deferribacterota", sex == "M", timepoint == "w5", diet == "Control Vehicle")

ps_pup_phylum_z_mean <- ps_pup_phylum_z %>%
  select(Phylum, z_score) %>%
  group_by(Phylum) %>%
  summarize(mean_z = mean(z_score, na.rm = TRUE))

ps_pup_phylum_check <- ps_pup_phylum %>%
  select(Sample, diet, sex, timepoint, Phylum, Abundance) %>%
  filter(Phylum == "Bacteroidota")

ps_pup_phylum %>%
  filter(Phylum == "Firmicutes") %>%
  ggplot(aes(x = Abundance)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  ggtitle("Distribution of Firmicutes Abundance") +
  theme_bw()

#compare mean and median
firmicutes_stats <- ps_pup_phylum %>%
  filter(Phylum == "Firmicutes") %>%
  summarise(mean_abundance = mean(Abundance), median_abundance = median(Abundance))
print(firmicutes_stats)

ps_pup_phylum_z %>%
  filter(Phylum == "Firmicutes") %>%
  ggplot(aes(x = z_score)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  ggtitle("Distribution of Firmicutes Z-Scores") +
  theme_bw()

mean(ps_pup_phylum_z$Abundance, na.rm = TRUE)
sd(ps_pup_phylum_z$Abundance, na.rm = TRUE)

#ps_pup_phylum_z_check <- ps_pup_phylum_z %>%
  #select(Sample, diet, sex, timepoint, Phylum, Abundance, z_score, mean, SD) %>%
  #filter(Phylum == "Cyanobacteria") #diet == "HFHS B. longum APC1472")

range(ps_pup_phylum_z$z_score, na.rm = TRUE)

###heatmap
##get rid of cyanobacteria because its all grey.
ps_pup_phylum_z <- ps_pup_phylum_z %>%
  filter(Phylum != "Cyanobacteria")

p.heat1 <- ggplot(ps_pup_phylum_z, aes(x = diet, y = Phylum)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Phylum") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2.2, 1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank(),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  facet_grid(sex ~ timepoint) + 
  ylab("Phylum") + 
  theme(panel.grid=element_blank()) +
  geom_hline(yintercept = seq(0.5, length(unique(ps_pup_phylum_z$Phylum)) - 0.5), color = "white")
p.heat1
getwd()
setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance")
ggsave("heatmap_group_byphylum_pups.png", p.heat1, width = 10, height = 6, dpi = 1200)
ggsave("heatmap_group_byphylum_pups.svg", p.heat1, width = 10, height = 6, device = "svg")

####Order
ps_pup_order <- PSr_pup %>%
  tax_glom(taxrank = "Order") %>%                     
  psmelt() %>%                                         
  #filter(Abundance > 0.02) %>%          
  arrange(Order)

#specific_orders <- c("Staphylococcales", "Eubacteriales", "Enterobacterales")

top_orders <- ps_pup_order %>%
  group_by(Order) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:13) %>%
  pull(Order)

#top_orders <- unique(c(top_orders, specific_orders))

ps_pup_order <- ps_pup_order %>%
  mutate(Order = ifelse(Order %in% top_orders, Order, "Other")) %>%
  mutate(Order = factor(Order, levels = c(top_orders, "Other")))

ps_pup_order_check <- ps_pup_order %>%
  select(Sample, diet, sex, timepoint, Order, Abundance) %>%
  filter(Order == "Enterobacterales")

###heatmap
ps_pup_order_z <- ps_pup_order %>%
  group_by(Order, sex, timepoint) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

range(ps_pup_order_z$z_score, na.rm = TRUE)

p.heat2 <- ggplot(ps_pup_order_z, aes(x = diet, y = Order)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Order") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2.2, 1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank(),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  facet_grid(sex ~ timepoint) + 
  ylab("Order") +
  geom_hline(yintercept = seq(0.5, length(unique(ps_pup_order_z$Order)) - 0.5), color = "white") +
  theme(panel.grid=element_blank())
p.heat2
ggsave("heatmap_group_byorder_pups.png", p.heat2, width = 10, height = 6, dpi = 1200)
ggsave("heatmap_group_byorder_pups.svg", p.heat2, width = 10, height = 6, device = "svg")

#####Family
###bars
ps_pup_family <- PSr_pup %>%
  tax_glom(taxrank = "Family") %>%                    
  psmelt() %>%                                        
  #filter(Abundance > 0.02) %>%           
  arrange(Family)

ps_pup_family_check <- ps_pup_family %>%
  filter(Family %in% c("Butyricicoccaceae", "Clostridiaceae", "Enterobacteriaceae",
                       "Enterococcaceae", "Streptococcaceae"))

#specific_families <- c("Butyricicoccaceae", "Clostridiaceae", "Enterobacteriaceae",
 #                      "Enterococcaceae", "Streptococcaceae")

specific_families <- "Butyricicoccaceae"

top_families <- ps_pup_family %>%
  group_by(Family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:14) %>%
  pull(Family)

top_families <- unique(c(top_families, specific_families))

ps_pup_family <- ps_pup_family %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other")) %>%  
  mutate(Family = factor(Family, levels = c(top_families, "Other")))  

###heatmap
ps_pup_family_z <- ps_pup_family %>%
  group_by(Family, sex, timepoint) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat3 <- ggplot(ps_pup_family_z, aes(x = diet, y = Family)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Family") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank(),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  facet_grid(sex ~ timepoint) + 
  ylab("Family") +
  geom_hline(yintercept = seq(0.5, length(unique(ps_pup_family_z$Family)) - 0.5), color = "white") +
  theme(panel.grid=element_blank())
p.heat3
ggsave("heatmap_group_byfamily_pups.png", p.heat3, width = 10, height = 6, dpi = 1200)
ggsave("heatmap_group_byfamily_pups.svg", p.heat3, width = 10, height = 6, device = "svg")

#####Genus
###bars
ps_pup_genus <- PSr_pup %>%
  tax_glom(taxrank = "Genus") %>%                     
  psmelt() %>%                                         
  #filter(Abundance > 0.02) %>%                        
  arrange(Genus) %>%
  mutate(Genus = gsub("Lachnospiraceae NK4A136 group", "Lachnospiraceae NK4A136", Genus),  # Shorten taxa names
         Genus = gsub("Rikenellaceae RC9 gut group", "Rikenellaceae RC9", Genus))  # Shorten taxa names

ps_pup_genus_check <- ps_pup_genus %>%
  filter(Genus %in% c("Acetatifactor", "Butyricicoccus", "Enterococcus", "Parabacteroides", "Roseburia"))

specific_genera <- c("Acetatifactor", "Butyricicoccus", "Enterococcus", "Parabacteroides", "Roseburia")

#specific_genera <- c("Acetatifactor", "Butyricicoccus", "Parabacteroides", "Roseburia")

top_genera <- ps_pup_genus %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:10) %>%
  pull(Genus)

top_genera <- unique(c(top_genera, specific_genera))

ps_pup_genus <- ps_pup_genus %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  mutate(Genus = factor(Genus, levels = c(top_genera, "Other")))

###heatmap
ps_pup_genus_z <- ps_pup_genus %>%
  group_by(Genus, sex, timepoint) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat4 <- ggplot(ps_pup_genus_z, aes(x = diet, y = Genus)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Genus") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank(),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  facet_grid(sex ~ timepoint) + 
  ylab("Genus") +
  geom_hline(yintercept = seq(0.5, length(unique(ps_pup_genus_z$Genus)) - 0.5), color = "white") +
  theme(panel.grid=element_blank())
p.heat4
ggsave("heatmap_group_bygenus_pups.png", p.heat4, width = 10, height = 6, dpi = 1200)
ggsave("heatmap_group_bygenus_pups.svg", p.heat4, width = 10, height = 6, device = "svg")



######################dams
#####Phylum
###bars
ps_dam_phylum <- PSr_dam %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level 
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum) 

ps_dam_phylum$Phylum <- factor(ps_dam_phylum$Phylum, levels = ps_dam_phylum %>%
                                 group_by(Phylum) %>%
                                 summarise(total_abundance = sum(Abundance)) %>%
                                 arrange(desc(total_abundance)) %>%
                                 pull(Phylum))

###Heatmap
ps_dam_phylum_z <- ps_dam_phylum %>%
  group_by(Phylum) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

#get rid of cyanobacteria
ps_dam_phylum_z <- ps_dam_phylum_z %>%
  filter(Phylum != "Cyanobacteria")

p.heat6 <- ggplot(ps_dam_phylum_z, aes(x = diet, y = Phylum)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Phylum") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2.2, 1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank(),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  #facet_grid(sex ~ timepoint) + 
  ylab("Phylum") +
  geom_hline(yintercept = seq(0.5, length(unique(ps_dam_phylum_z$Phylum)) - 0.5), color = "white") +
  theme(panel.grid=element_blank())
p.heat6
ggsave("heatmap_group_byphylum_dams.png", p.heat6, width = 10, height = 6, dpi = 1200)
ggsave("heatmap_group_byphylum_dams.svg", p.heat6, width = 10, height = 6, device = "svg")

####Order
###bars
ps_dam_order <- PSr_dam %>%
  tax_glom(taxrank = "Order") %>%                     
  psmelt() %>%                                         
  #filter(Abundance > 0.02) %>%                        
  arrange(Order)

ps_dam_order_check <- ps_dam_order %>%
  filter(Order %in% c("Staphylococcales", "Eubacteriales", "Enterobacterales"))

#specific_orders <- c("Staphylococcales", "Eubacteriales", "Enterobacterales")

specific_orders <- "Clostridiales"

top_orders <- ps_dam_order %>%
  group_by(Order) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:12) %>%
  pull(Order)

top_orders <- unique(c(top_orders, specific_orders))

ps_dam_order <- ps_dam_order %>%
  mutate(Order = ifelse(Order %in% top_orders, Order, "Other")) %>%
  mutate(Order = factor(Order, levels = c(top_orders, "Other")))

###heatmap
ps_dam_order_z <- ps_dam_order %>%
  group_by(Order) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat7 <- ggplot(ps_dam_order_z, aes(x = diet, y = Order)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Order") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2.2, 1)) +
  theme_bw() +
  #theme(panel.grid=element_blank()) +
  theme(
    axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank(),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  #facet_grid(sex ~ timepoint) + 
  ylab("Order") +
  geom_hline(yintercept = seq(0.5, length(unique(ps_dam_order_z$Order)) - 0.5), color = "white") +
  theme(panel.grid=element_blank())
p.heat7
ggsave("heatmap_group_byorder_dams.png", p.heat7, width = 10, height = 6, dpi = 1200)
ggsave("heatmap_group_byorder_dams.svg", p.heat7, width = 10, height = 6, device = "svg")

#####Family
###bars
ps_dam_family <- PSr_dam %>%
  tax_glom(taxrank = "Family") %>%                    
  psmelt() %>%                                        
  #filter(Abundance > 0.02) %>%                     
  arrange(Family)

ps_dam_family_check <- ps_dam_family %>%
  filter(Family %in% c("Butyricicoccaceae", "Clostridiaceae", "Enterobacteriaceae",
                       "Enterococcaceae", "Streptococcaceae"))

#specific_families <- c("Butyricicoccaceae", "Clostridiaceae", "Enterobacteriaceae",
 #                     "Enterococcaceae", "Streptococcaceae")

specific_families <- "Butyricicoccaceae"

top_families <- ps_dam_family %>%
  group_by(Family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:14) %>%
  pull(Family)

top_families <- unique(c(top_families, specific_families))

ps_dam_family <- ps_dam_family %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other")) %>%  
  mutate(Family = factor(Family, levels = c(top_families, "Other")))  

###heatmap
ps_dam_family_z <- ps_dam_family %>%
  group_by(Family) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat8 <- ggplot(ps_dam_family_z, aes(x = diet, y = Family)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Family") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank(),
    strip.background = element_rect(colour="black", fill="white")
  ) +
 #facet_grid(sex ~ timepoint) + 
  ylab("Family") +
  geom_hline(yintercept = seq(0.5, length(unique(ps_dam_family_z$Family)) - 0.5), color = "white") +
  theme(panel.grid=element_blank())
p.heat8
ggsave("heatmap_group_byfamily_dams.png", p.heat8, width = 10, height = 6, dpi = 1200)
ggsave("heatmap_group_byfamily_dams.svg", p.heat8, width = 10, height = 6, device = "svg")

#####Genus
###bars
ps_dam_genus <- PSr_dam %>%
  tax_glom(taxrank = "Genus") %>%                     
  psmelt() %>%                                         
  #filter(Abundance > 0.02) %>%                                  
  arrange(Genus) %>%
  mutate(Genus = gsub("Lachnospiraceae NK4A136 group", "Lachnospiraceae NK4A136", Genus),  # Shorten taxa names
         Genus = gsub("Rikenellaceae RC9 gut group", "Rikenellaceae RC9", Genus))  # Shorten taxa names

ps_dam_genus_check <- ps_dam_genus %>%
  filter(Genus %in% c("Acetatifactor", "Butyricicoccus", "Enterococcus", "Parabacteroides", "Roseburia"))

specific_genera <- c("Acetatifactor", "Butyricicoccus", "Enterococcus", "Parabacteroides", "Roseburia")

top_genera <- ps_dam_genus %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:10) %>%
  pull(Genus)

top_genera <- unique(c(top_genera, specific_genera))

ps_dam_genus <- ps_dam_genus %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  mutate(Genus = factor(Genus, levels = c(top_genera, "Other")))

###heatmap
ps_dam_genus_z <- ps_dam_genus %>%
  group_by(Genus) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat9 <- ggplot(ps_dam_genus_z, aes(x = diet, y = Genus)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Genus") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank(),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  #facet_grid(sex ~ timepoint) + 
  ylab("Genus") +
  geom_hline(yintercept = seq(0.5, length(unique(ps_dam_genus_z$Genus)) - 0.5), color = "white") +
  theme(panel.grid=element_blank())
p.heat9
ggsave("heatmap_group_bygenus_dams.png", p.heat9, width = 10, height = 6, dpi = 1200)
ggsave("heatmap_group_bygenus_dams.svg", p.heat9, width = 10, height = 6, device = "svg")
