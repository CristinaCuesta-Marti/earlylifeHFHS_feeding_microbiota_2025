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
###Bar charts
ps_pup_phylum <- PSr_pup %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level 
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum) 
str(ps_pup_phylum)

ps_pup_phylum$Phylum <- factor(ps_pup_phylum$Phylum, levels = ps_pup_phylum %>%
                              group_by(Phylum) %>%
                              summarise(total_abundance = sum(Abundance)) %>%
                              arrange(desc(total_abundance)) %>%
                              pull(Phylum))

plot1 <- ggplot(ps_pup_phylum, aes(x = diet, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance pups - Phylum") +
  facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot1


###Heatmap

ps_pup_phylum_z <- ps_pup_phylum %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

p.heat1 <- ggplot(ps_pup_phylum_z, aes(x = diet, y = Phylum)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Phylum") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 2)) +
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
  geom_hline(yintercept = seq(0.5, length(unique(ps_pup_phylum_z$Phylum)) - 0.5), color = "white")

####Order
###bars
ps_pup_order <- PSr_pup %>%
  tax_glom(taxrank = "Order") %>%                     
  psmelt() %>%                                         
  filter(Abundance > 0.02) %>%                         
  arrange(Order)

top_orders <- ps_pup_order %>%
  group_by(Order) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:12) %>%
  pull(Order)

ps_pup_order <- ps_pup_order %>%
  mutate(Order = ifelse(Order %in% top_orders, Order, "Other")) %>%
  mutate(Order = factor(Order, levels = c(top_orders, "Other")))

plot2 <- ggplot(ps_pup_order, aes(x = diet, y = Abundance, fill = Order)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance pups - Order") +
  facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot2

###heatmap
ps_pup_order_z <- ps_pup_order %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

p.heat2 <- ggplot(ps_pup_order_z, aes(x = diet, y = Order)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Order") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 2)) +
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

#####Family
###bars
ps_pup_family <- PSr_pup %>%
  tax_glom(taxrank = "Family") %>%                    
  psmelt() %>%                                        
  filter(Abundance > 0.02) %>%                         
  arrange(Family)

top_families <- ps_pup_family %>%
  group_by(Family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:12) %>%
  pull(Family)

ps_pup_family <- ps_pup_family %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other")) %>%  
  mutate(Family = factor(Family, levels = c(top_families, "Other")))  

plot3 <- ggplot(ps_pup_family, aes(x = diet, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance pups - Family") +
  facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

###heatmap
ps_pup_family_z <- ps_pup_family %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

p.heat3 <- ggplot(ps_pup_family_z, aes(x = diet, y = Family)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Family") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 2)) +
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

#####Genus
###bars
ps_pup_genus <- PSr_pup %>%
  tax_glom(taxrank = "Genus") %>%                     
  psmelt() %>%                                         
  filter(Abundance > 0.02) %>%                         
  arrange(Genus) %>%
  mutate(Genus = gsub("Lachnospiraceae NK4A136 group", "Lachnospiraceae NK4A136", Genus),  # Shorten taxa names
         Genus = gsub("Rikenellaceae RC9 gut group", "Rikenellaceae RC9", Genus))  # Shorten taxa names

top_genera <- ps_pup_genus %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:15) %>%
  pull(Genus)

ps_pup_genus <- ps_pup_genus %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  mutate(Genus = factor(Genus, levels = c(top_genera, "Other")))

plot4 <- ggplot(ps_pup_genus, aes(x = diet, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance pups - Genus") +
  facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot4

###heatmap
ps_pup_genus_z <- ps_pup_genus %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

p.heat4 <- ggplot(ps_pup_genus_z, aes(x = diet, y = Genus)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Genus") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 2)) +
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

#####species
##Pups
view(tax_table)

ps_pup_species <- PSr_pup %>%
  tax_glom(taxrank = "Species") %>%                     
  psmelt() %>%    
  mutate(Genus_species = paste(Genus, Species, sep = "_")) %>%
  #filter(Abundance > 0.02) %>%                        
  arrange(Genus_species) 
  #mutate(Genus = gsub("Lachnospiraceae NK4A136 group", "Lachnospiraceae NK4A136", Genus),  # Shorten taxa names
   #      Genus = gsub("Rikenellaceae RC9 gut group", "Rikenellaceae RC9", Genus))  # Shorten taxa names

top_species <- ps_pup_species %>%
  group_by(Genus_species) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:15) %>%
  pull(Genus_species)

ps_pup_species <- ps_pup_species %>%
  mutate(Genus_species = ifelse(Genus_species %in% top_species, Genus_species, "Other")) %>%
  mutate(Genus_species = factor(Genus_species, levels = c(top_species, "Other")))

plot_pup_species <- ggplot(ps_pup_species, aes(x = diet, y = Abundance, fill = Genus_species)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance pups - Species") +
  facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot_pup_species

#### B. longum
diet_colours <- c("Control Vehicle" = "#a0a0a4", "HFHS Vehicle" = "#f94040", "HFHS FOS+GOS" = "#addead", "HFHS B. longum APC1472" = "#e1c180")

ps_pup_longum <- PSr_pup %>%
  tax_glom(taxrank = "Species") %>%                     
  psmelt() %>%  
  mutate(Genus_species = paste(Genus, Species, sep = "_")) %>%
  #filter(Abundance > 0.02) %>%                         
  arrange(Species)  
#Select species of interest
species_of_interest <- "Bifidobacterium_longum"
# Filter the data to include only the species of interest
filtered_data <- ps_pup_longum %>%
  filter(Genus_species == species_of_interest)
# Define the color for the species of interest
species_colors <- setNames("#DA5724", species_of_interest)
# Create a named vector for the labels with line breaks
species_labels <- setNames(str_wrap(species_of_interest, width = 20), species_of_interest)

plot5 <- ggplot(filtered_data, aes(x = diet, y = Abundance, fill = diet)) + 
  geom_bar(stat = "identity") + 
  labs(x="", y="Relative abundance", title = expression("Relative Abundance pups - " * italic("B. longum"))) +
  facet_grid(sex ~ timepoint, scales= "free_x") +
  scale_fill_manual(values = diet_colours) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot5

#heatmap_longum
filtered_data_z <- filtered_data %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

plot5_z <- ggplot(filtered_data_z, aes(x = diet, y = z_score, fill = diet)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = diet_colours) +
  labs(
    x = "", 
    y = "z_score", 
    title = expression("z-scores rel. abund. pups - " * italic("B. longum"))
  ) +
  facet_grid(sex ~ timepoint, scales = "free_x") +
  theme_bw() + 
  scale_y_continuous(limits = c(-2, 2)) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.text = element_text(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
plot5_z

plot5_z_nolimits <- ggplot(filtered_data_z, aes(x = diet, y = z_score, fill = diet)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = diet_colours) +
  labs(
    x = "", 
    y = "z_score", 
    title = expression("z-scores rel. abund. pups - " * italic("B. longum"))
  ) +
  facet_grid(sex ~ timepoint, scales = "free_x") +
  theme_bw() + 
  scale_y_continuous() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.text = element_text(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
plot5_z_nolimits

plot5_z_limits1025 <- ggplot(filtered_data_z, aes(x = diet, y = z_score, fill = diet)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = diet_colours) +
  labs(
    x = "", 
    y = "z_score", 
    title = expression("z-scores rel. abund. pups - " * italic("B. longum"))
  ) +
  facet_grid(sex ~ timepoint, scales = "free_x") +
  theme_bw() + 
  scale_y_continuous(limits = c(-10, 25)) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.text = element_text(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
plot5_z_limits1025

plot5_z_limits310 <- ggplot(filtered_data_z, aes(x = diet, y = z_score, fill = diet)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = diet_colours) +
  labs(
    x = "", 
    y = "z_score", 
    title = expression("z-scores rel. abund. pups - " * italic("B. longum"))
  ) +
  facet_grid(sex ~ timepoint, scales = "free_x") +
  theme_bw() + 
  scale_y_continuous(limits = c(-3, 10)) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.text = element_text(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
plot5_z_limits310

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

plot6 <- ggplot(ps_dam_phylum, aes(x = diet, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance dams - Phylum") +
  #facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot6

###Heatmap
ps_dam_phylum_z <- ps_dam_phylum %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

p.heat6 <- ggplot(ps_dam_phylum_z, aes(x = diet, y = Phylum)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Phylum") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 2)) +
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
  geom_hline(yintercept = seq(0.5, length(unique(ps_dam_phylum_z$Phylum)) - 0.5), color = "white")
p.heat6

####Order
###bars
ps_dam_order <- PSr_dam %>%
  tax_glom(taxrank = "Order") %>%                     
  psmelt() %>%                                         
  filter(Abundance > 0.02) %>%                         
  arrange(Order)

top_orders <- ps_dam_order %>%
  group_by(Order) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:12) %>%
  pull(Order)

ps_dam_order <- ps_dam_order %>%
  mutate(Order = ifelse(Order %in% top_orders, Order, "Other")) %>%
  mutate(Order = factor(Order, levels = c(top_orders, "Other")))

plot7 <- ggplot(ps_dam_order, aes(x = diet, y = Abundance, fill = Order)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance dams - Order") +
  #facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot7

###heatmap
ps_dam_order_z <- ps_dam_order %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

p.heat7 <- ggplot(ps_dam_order_z, aes(x = diet, y = Order)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Order") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 2)) +
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

#####Family
###bars
ps_dam_family <- PSr_dam %>%
  tax_glom(taxrank = "Family") %>%                    
  psmelt() %>%                                        
  filter(Abundance > 0.02) %>%                         
  arrange(Family)

top_families <- ps_dam_family %>%
  group_by(Family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:12) %>%
  pull(Family)

ps_dam_family <- ps_dam_family %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other")) %>%  
  mutate(Family = factor(Family, levels = c(top_families, "Other")))  

plot8 <- ggplot(ps_dam_family, aes(x = diet, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance dams - Family") +
  #facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot8

###heatmap
ps_dam_family_z <- ps_dam_family %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

p.heat8 <- ggplot(ps_dam_family_z, aes(x = diet, y = Family)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Family") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 2)) +
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

#####Genus
###bars
ps_dam_genus <- PSr_dam %>%
  tax_glom(taxrank = "Genus") %>%                     
  psmelt() %>%                                         
  filter(Abundance > 0.02) %>%                         
  arrange(Genus) %>%
  mutate(Genus = gsub("Lachnospiraceae NK4A136 group", "Lachnospiraceae NK4A136", Genus),  # Shorten taxa names
         Genus = gsub("Rikenellaceae RC9 gut group", "Rikenellaceae RC9", Genus))  # Shorten taxa names

top_genera <- ps_dam_genus %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:15) %>%
  pull(Genus)

ps_dam_genus <- ps_dam_genus %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  mutate(Genus = factor(Genus, levels = c(top_genera, "Other")))

plot9 <- ggplot(ps_dam_genus, aes(x = diet, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance dams - Genus") +
  #facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot9

###heatmap
ps_dam_genus_z <- ps_dam_genus %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

p.heat9 <- ggplot(ps_dam_genus_z, aes(x = diet, y = Genus)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Genus") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 2)) +
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

####Species
##Dams
ps_dam_species <- PSr_dam %>%
  tax_glom(taxrank = "Species") %>%                     
  psmelt() %>%    
  mutate(Genus_species = paste(Genus, Species, sep = "_")) %>%
  #filter(Abundance > 0.02) %>%                        
  arrange(Genus_species) 
#mutate(Genus = gsub("Lachnospiraceae NK4A136 group", "Lachnospiraceae NK4A136", Genus),  # Shorten taxa names
#      Genus = gsub("Rikenellaceae RC9 gut group", "Rikenellaceae RC9", Genus))  # Shorten taxa names

top_species <- ps_dam_species %>%
  group_by(Genus_species) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:15) %>%
  pull(Genus_species)

ps_dam_species <- ps_dam_species %>%
  mutate(Genus_species = ifelse(Genus_species %in% top_species, Genus_species, "Other")) %>%
  mutate(Genus_species = factor(Genus_species, levels = c(top_species, "Other")))

plot_dam_species <- ggplot(ps_dam_species, aes(x = diet, y = Abundance, fill = Genus_species)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance dams - Species") +
  #facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot_dam_species

#### B. longum
ps_dam_longum <- PSr_dam %>%
  tax_glom(taxrank = "Species") %>%                     
  psmelt() %>%  
  mutate(Genus_species = paste(Genus, Species, sep = "_")) %>%
  #filter(Abundance > 0.02) %>%                         
  arrange(Species)  
#Select species of interest
species_of_interest <- "Bifidobacterium_longum"
# Filter the data to include only the species of interest
filtered_datad <- ps_dam_longum %>%
  filter(Genus_species == species_of_interest)
# Define the color for the species of interest
species_colors <- setNames("#DA5724", species_of_interest)
# Create a named vector for the labels with line breaks
species_labels <- setNames(str_wrap(species_of_interest, width = 20), species_of_interest)

plot10 <- ggplot(filtered_datad, aes(x = diet, y = Abundance, fill = diet)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = diet_colours) +
  labs(x="", y="Relative abundance", title = expression("Relative Abundance dams - " * italic("B. longum"))) +
  #facet_grid(sex ~ timepoint, scales= "free_x") +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot10

##Dams longum z-score
filtered_datad_z <- filtered_datad %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

plot10_z <- ggplot(filtered_datad_z, aes(x = diet, y = z_score, fill = diet)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = diet_colours) +
  labs(
    x = "", 
    y = "z_score", 
    title = expression("z-scores rel. abund. dams - " * italic("B. longum"))
  ) +
  #facet_grid(sex ~ timepoint, scales = "free_x") +
  theme_bw() + 
  scale_y_continuous(limits = c(-2, 2)) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.text = element_text(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
plot10_z

plot10_z_nolimits <- ggplot(filtered_datad_z, aes(x = diet, y = z_score, fill = diet)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = diet_colours) +
  labs(
    x = "", 
    y = "z_score", 
    title = expression("z-scores rel. abund. dams - " * italic("B. longum"))
  ) +
  #facet_grid(sex ~ timepoint, scales = "free_x") +
  theme_bw() + 
  scale_y_continuous() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.text = element_text(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
plot10_z_nolimits

plot10_z_limits1025 <- ggplot(filtered_datad_z, aes(x = diet, y = z_score, fill = diet)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = diet_colours) +
  labs(
    x = "", 
    y = "z_score", 
    title = expression("z-scores rel. abund. dams - " * italic("B. longum"))
  ) +
  #facet_grid(sex ~ timepoint, scales = "free_x") +
  theme_bw() + 
  scale_y_continuous(limits = c(-10, 25)) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.text = element_text(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
plot10_z_limits1025

plot10_z_limits310 <- ggplot(filtered_datad_z, aes(x = diet, y = z_score, fill = diet)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = diet_colours) +
  labs(
    x = "", 
    y = "z_score", 
    title = expression("z-scores rel. abund. dams - " * italic("B. longum"))
  ) +
  #facet_grid(sex ~ timepoint, scales = "free_x") +
  theme_bw() + 
  scale_y_continuous(limits = c(-3, 10)) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.text = element_text(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
plot10_z_limits310


######same with ggpubr. I think this way is better.No subsetting. Pups
stat_pup_phylum <- compare_means(data = ps_pup_phylum, 
                                 Abundance ~ diet, 
                                 group.by = c("sex", "timepoint", "Phylum"))

stat_pup_order <- compare_means(data = ps_pup_order, 
                                Abundance ~ diet, 
                                group.by = c("sex", "timepoint", "Order"))

stat_pup_family <- compare_means(data = ps_pup_family, 
                                 Abundance ~ diet, 
                                 group.by = c("sex", "timepoint", "Family"))

stat_pup_genus <- compare_means(data = ps_pup_genus, 
                                Abundance ~ diet, 
                                group.by = c("sex", "timepoint", "Genus"))

stat_pup_species <- compare_means(data = ps_pup_species, 
                                Abundance ~ diet, 
                                group.by = c("sex", "timepoint", "Genus_species"))  

stat_pup_longum <- compare_means(data = filtered_data, 
                                 Abundance ~ diet, 
                                 group.by = c("sex", "timepoint", "Genus_species"))

write_xlsx(stat_pup_phylum, "stats_relativeabundance_pups_phylum.xlsx")
write_xlsx(stat_pup_order, "stats_relativeabundance_pups_order.xlsx")
write_xlsx(stat_pup_family, "stats_relativeabundance_pups_family.xlsx")
write_xlsx(stat_pup_genus, "stats_relativeabundance_pups_genus.xlsx")
write_xlsx(stat_pup_longum, "stats_relativeabundance_pups_longum.xlsx")

######same with ggpubr. I think this way is better.No subsetting. Dams
stat_dam_phylum <- compare_means(data = ps_dam_phylum, 
                                 Abundance ~ diet,
                                 group.by = "Phylum")

stat_dam_order <- compare_means(data = ps_dam_order, 
                                Abundance ~ diet,
                                group.by = "Order")

stat_dam_family <- compare_means(data = ps_dam_family, 
                                 Abundance ~ diet,
                                 group.by = "Family")

stat_dam_genus <- compare_means(data = ps_dam_genus, 
                                Abundance ~ diet,
                                group.by = "Genus")

stat_dam_species <- compare_means(data = ps_dam_species, 
                                Abundance ~ diet, 
                                group.by = c("Genus_species")) 

stat_dam_longum <- compare_means(data = filtered_datad, 
                                 Abundance ~ diet,
                                 group.by = "Genus_species")

write_xlsx(stat_dam_phylum, "stats_relativeabundance_dams_phylum.xlsx")
write_xlsx(stat_dam_family, "stats_relativeabundance_dams_family.xlsx")
write_xlsx(stat_dam_order, "stats_relativeabundance_dams_order.xlsx")
write_xlsx(stat_dam_genus, "stats_relativeabundance_dams_genus.xlsx")
write_xlsx(stat_dam_longum, "stats_relativeabundance_dams_longum.xlsx")

############Generate abundance tables for each taxon.
########Pups
phylum_abundance_table_pup <- ps_pup_phylum %>%
  group_by(Phylum, sex, timepoint, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

order_abundance_table_pup <- ps_pup_order %>%
  group_by(Order, sex, timepoint, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

family_abundance_table_pup <- ps_pup_family %>%
  group_by(Family, sex, timepoint, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

genus_abundance_table_pup <- ps_pup_genus %>%
  group_by(Genus, sex, timepoint, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

species_abundance_table_pup <- ps_pup_species %>%
  group_by(Genus_species, sex, timepoint, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

longum_abundance_table_pup <- filtered_data %>%
  group_by(Genus_species, sex, timepoint, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance")
write_xlsx(phylum_abundance_table_pup, "phylum_relativeabundance_pup.xlsx")
write_xlsx(order_abundance_table_pup, "order_relativeabundance_pup.xlsx")
write_xlsx(family_abundance_table_pup, "family_relativeabundance_pup.xlsx")
write_xlsx(genus_abundance_table_pup, "genus_relativeabundance_pup.xlsx")
write_xlsx(longum_abundance_table_pup, "longum_relativeabundance_pup.xlsx")

########Dams
phylum_abundance_table_dam <- ps_dam_phylum %>%
  group_by(Phylum, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

order_abundance_table_dam <- ps_dam_order %>%
  group_by(Order, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

family_abundance_table_dam <- ps_dam_family %>%
  group_by(Family, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

genus_abundance_table_dam <- ps_dam_genus %>%
  group_by(Genus, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

species_abundance_table_dam <- ps_dam_species %>%
  group_by(Genus_species, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

longum_abundance_table_dam <- filtered_datad %>%
  group_by(Genus_species, diet) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance")
write_xlsx(phylum_abundance_table_dam, "phylum_relativeabundance_dam.xlsx")
write_xlsx(order_abundance_table_dam, "order_relativeabundance_dam.xlsx")
write_xlsx(family_abundance_table_dam, "family_relativeabundance_dam.xlsx")
write_xlsx(genus_abundance_table_dam, "genus_relativeabundance_dam.xlsx")
write_xlsx(longum_abundance_table_dam, "longum_relativeabundance_dam.xlsx")