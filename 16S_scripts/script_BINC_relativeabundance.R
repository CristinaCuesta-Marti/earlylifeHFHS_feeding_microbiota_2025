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

setwd("C:/Users/Dara/OneDrive - University College Cork/ADATA onedrive/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items")

#Metadata or sample_data
samdf <- read_excel("Metadata_BINC_maxee24_ps.xlsx")
sam2 <- samdf %>% remove_rownames %>% column_to_rownames(var="sample")

#Split diet
sam2 <- sam2 %>%
  separate(diet, into = c("diet", "treatment"), sep = " ", extra = "merge") %>%
  mutate(#diet = factor(diet, levels = c("Control Vehicle", "HFHS Vehicle", "HFHS FOS+GOS", "HFHS B. longum APC1472")),
         diet = factor(diet, levels = c("Control", "HFHS")),
         treatment = factor(treatment, levels = c("Vehicle", "FOS+GOS", "B. longum APC1472")),
         #sex = factor(sex, levels = c("F", "M")),
         #timepoint = factor(timepoint, levels = c("w5", "w10")))
  )


sam2$sex <- factor(sam2$sex, levels = c("M", "F"), labels = c("Male", "Female"))
sam2$timepoint <- factor(sam2$timepoint, levels = c("w5", "w10"), labels = c("week 5", "week 10"))

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
#sample_data_reordered <- sample_data(ps3)
#sample_data_reordered <- as.data.frame(sample_data_reordered)
#sample_data_reordered$diet <- factor(sample_data_reordered$diet, levels = c("Control Vehicle", "HFHS Vehicle"
                                                                            #, "HFHS FOS+GOS", "HFHS B. longum APC1472"))
#sample_data_reordered$timepoint <- factor(sample_data_reordered$timepoint, levels = c("week 5", "week 10"))
#sample_data_reordered$sex <- factor(sample_data_reordered$sex, levels = c("Male", "Female"))
#sample_data <- sample_data_reordered
#ps3_reordered <- phyloseq(otu_table, tax_table, sample_data)
#sample_data(ps3_reordered)
#class(ps3_reordered)

##separate the pups and dams, and make do separate ps based on age column

ps_pup <- subset_samples(ps3, sample_data(ps)$age =="pup")
ps_dam <- subset_samples(ps3, sample_data(ps)$age =="dam")


#Make relative abundances
PSr_pup <- transform_sample_counts(ps_pup, function(x) x/sum(x))
PSr_dam <- transform_sample_counts(ps_dam, function(x) x/sum(x))

sample_data(PSr_pup)

sample_sums <- sample_sums(PSr_pup)
all(sample_sums == 1)

#########################pups
#####Phylum
###Bars
ps_pup_phylum <- PSr_pup %>%
  tax_glom(taxrank = "Phylum") %>%                    
  psmelt() %>%                                         
  filter(Abundance > 0.02) %>%  # turn this off for getting table, doing stats etc.             
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
getwd()
setwd("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance")
#ggsave("phylum_pups.png", plot1, dpi = 1200)
#ggsave("phylum_pups.svg", plot1, device = "svg")

###Heatmap
ps_pup_phylum_z <- ps_pup_phylum %>%
  group_by(Phylum) %>%
  mutate(mean = mean(Abundance)) %>%
  mutate(SD = sd(Abundance)) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance))

ps_pup_phylum_z_mean <- ps_pup_phylum_z %>%
  select(Phylum, z_score) %>%
  group_by(Phylum) %>%
  summarize(mean_z = mean(z_score, na.rm = TRUE))

#getwd()
#write_xlsx(ps_pup_phylum_z, "check.xlsx")

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

ps_pup_phylum_z_check <- ps_pup_phylum_z %>%
  select(Sample, diet, sex, timepoint, Phylum, Abundance, z_score, mean, SD) %>%
  filter(Phylum == "Deferribacterota", sex == "F", timepoint == "w10", diet == "HFHS Vehicle")

range(ps_pup_phylum_z$z_score, na.rm = TRUE)

p.heat1 <- ggplot(ps_pup_phylum_z, aes(x = diet, y = Phylum)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Phylum") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2.1, 1.2)) +
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
ggsave("heatmap_group_byphylum_pups.png", p.heat1, width = 10, height = 6, dpi = 1200)
ggsave("heatmap_group_byphylum_pups.svg", p.heat1, width = 10, height = 6, device = "svg")

####Order
###bars
ps_pup_order <- PSr_pup %>%
  tax_glom(taxrank = "Order") %>%                     
  psmelt() %>%                                         
  #filter(Abundance > 0.02) %>%          
  arrange(Order)

ps_pup_order_check <- ps_pup_order %>%
  filter(Order %in% c("Staphylococcales", "Eubacteriales", "Enterobacterales"))

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
ggsave("order_pups.png", plot2, dpi = 1200)
ggsave("order_pups.svg", plot2, device = "svg")

###heatmap
ps_pup_order_z <- ps_pup_order %>%
  group_by(Order) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

range(ps_pup_order_z$z_score, na.rm = TRUE)

p.heat2 <- ggplot(ps_pup_order_z, aes(x = diet, y = Order)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Order") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 0.8)) +
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
plot3
ggsave("family_pups.png", plot3, dpi = 1200)
ggsave("family_pups.svg", plot3, device = "svg")

###heatmap
ps_pup_family_z <- ps_pup_family %>%
  group_by(Family) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat3 <- ggplot(ps_pup_family_z, aes(x = diet, y = Family)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Family") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-2, 0.8)) +
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
ggsave("genus_pups.png", plot4, dpi = 1200)
ggsave("genus_pups.svg", plot4, device = "svg")

###heatmap
ps_pup_genus_z <- ps_pup_genus %>%
  group_by(Genus) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat4 <- ggplot(ps_pup_genus_z, aes(x = diet, y = Genus)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance pups - Genus") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-1.5, 0.5)) +
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
#getwd()
setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/Bar_charts")
ggsave("species_pups.png", plot_pup_species, dpi = 1200)
ggsave("species_pups.svg", plot_pup_species, device = "svg")

#### B. longum #this plot is not actually correct
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
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot5
#getwd()
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance")
ggsave("longum_pups.png", plot5, dpi = 1200)
ggsave("longum_pups.svg", plot5, device = "svg")

#heatmap_longum
filtered_data_z <- filtered_data %>%
  group_by(Genus_species) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

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
ggsave("longum_pups_zscore.png", plot5_z, dpi = 1200)
ggsave("longum_pups_zscore.svg", plot5_z, device = "svg")

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
ggsave("longum_pups_zscore_nolimits.png", plot5_z_nolimits, dpi = 1200)
ggsave("longum_pups_zscore_nolimits.svg", plot5_z_nolimits, device = "svg")

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
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/Bar_charts")
ggsave("longum_pups_zscore_limits1025.png", plot5_z_limits1025, dpi = 1200)
ggsave("longum_pups_zscore_limits1025.svg", plot5_z_limits1025, device = "svg")

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
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/Bar_charts")
ggsave("longum_pups_zscore_limits310.png", plot5_z_limits310, dpi = 1200)
ggsave("longum_pups_zscore_limits310.svg", plot5_z_limits310, device = "svg")

##########Bifidobacterium species bar chart. Bif species relative to one another.
ps_pup_species <- PSr_pup %>%
  tax_glom(taxrank = "Species") %>%                     
  psmelt() %>%    
  mutate(Genus_species = paste(Genus, Species, sep = "_")) %>%
  arrange(Genus_species) 

ps_pup_bif <- ps_pup_species %>%
  filter(Genus == "Bifidobacterium") 
getwd()
#write_xlsx(ps_pup_bif, "dataframe_pups_Bifidobacterium.xlsx")

#also send her absolute abundance for bif
ps_absolute_pup_species <- ps_pup %>%
  tax_glom(taxrank = "Species") %>%                     
  psmelt() %>%    
  mutate(Genus_species = paste(Genus, Species, sep = "_")) %>%
  #filter(Abundance > 0.02) %>%                        
  arrange(Genus_species) 

ps_absolute_pup_bif <- ps_absolute_pup_species %>%
  filter(Genus == "Bifidobacterium")
write_xlsx(ps_absolute_pup_bif, "dataframe_absolute_pups_Bifidobacterium.xlsx")

plot_pup_bif <- ggplot(ps_pup_species, aes(x = diet, y = Abundance, fill = ifelse(
  Genus == "Bifidobacterium", Genus_species, NA))) + 
  geom_bar(stat = "identity", position = "identity") + 
  labs(x="", y="Relative abundance", title = "Relative Abundance pups - Bifidobacterium") +
  facet_grid(sex ~ timepoint, scales= "free_x") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot_pup_bif


#getwd()
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/Bar_charts")
#ggsave("species_pups.png", plot_pup_species, dpi = 1200)
#ggsave("species_pups.svg", plot_pup_species, device = "svg")

#################Cristina's script for Bif species
diet_colours <- c("Control Vehicle" = "#a0a0a4", "HFHS Vehicle" = "#f94040", "HFHS FOS+GOS" = "#addead", "HFHS B. longum APC1472" = "#e1c180")

ps_pup_pseudolongum <- ps_pup_bif %>%
  filter(Genus_species == "Bifidobacterium_pseudolongum")
ps_pup_longum <- ps_pup_bif %>%
  filter(Genus_species == "Bifidobacterium_longum")
ps_pup_breve <- ps_pup_bif %>%
  filter(Genus_species == "Bifidobacterium_breve")

# comparing mean relative abundances
Bifo_mean <- ps_pup_bif %>%
  group_by(timepoint, sex, Genus_species, diet) %>%
  summarise(mean_abundance = mean(Abundance),
            se = sd(Abundance)/sqrt(n()),
            .groups = 'drop'
            )

plot_abundant_bif <- Bifo_mean %>%
  ggplot(aes(x = mean_abundance, y = Genus_species)) + 
  geom_errorbarh(aes(xmin = mean_abundance - se,
                     xmax = mean_abundance + se,
                     color = diet),
                 height = 0.2,
                 position = position_dodge(width = -0.7)) + # Changed to negative value
  # Add point at mean
  geom_point(aes(color = diet),
             size = 3,
             position = position_dodge(width = -0.7)) + # Changed to negative value
  scale_color_manual(values = diet_colours) +
  facet_grid(timepoint ~ sex) +
  xlim(0, 0.3) +
  theme_bw() +
  ylab(NULL) +
  xlab("Relative abundance") +
  theme(legend.position = "bottom",
        legend.background = element_rect(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
plot_abundant_bif
ggsave("pups_bifidobacterium_dots.png", plot_abundant_bif, width = 10, height = 6, dpi = 1200)
ggsave("pups_bifidobacterium_dots.svg", plot_abundant_bif, width = 10, height = 6, device = "svg")

# pseudolongum
plot_pup_pseudolongum <- ps_pup_pseudolongum %>%
  ggplot(aes(x = sex, y = Abundance, fill = diet, colour = diet)) +
  # Change to boxplot
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
              size = 2, colour = "black") +
  geom_vline(xintercept = 1.5, colour = "darkgray", linetype = 'dashed') +
  ylim(0, 0.34) +  # Add this line to set y-axis limits
  facet_wrap(~timepoint, scales = "fixed") +  # Changed to "fixed" to maintain same scale
  scale_fill_manual(values = diet_colours) +
  theme_bw() +
  labs(x="", y="Relative abundance", title = expression("Relative abundance of " * italic("B. pseudolongum"))) 
plot_pup_pseudolongum
#setwd("E:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/fancy_new_plots")
ggsave("pups_pseudolongum.png", plot_pup_pseudolongum, dpi = 1200)
ggsave("pups_pseudolongum.svg", plot_pup_pseudolongum, device = "svg")

# longum
plot_pup_longum <- ps_pup_longum %>%
  ggplot(aes(x = sex, y = Abundance, fill = diet, colour = diet)) +
  # Change to boxplot
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
              size = 2, colour = "black") +
  geom_vline(xintercept = 1.5, colour = "darkgray", linetype = 'dashed') +
  ylim(0, 0.02) +  # Add this line to set y-axis limits
  facet_wrap(~timepoint, scales = "fixed") +  # Changed to "fixed" to maintain same scale
  scale_fill_manual(values = diet_colours) +
  theme_bw() +
  labs(x="", y="Relative abundance", title = expression("Relative abundance of " * italic("B. longum"))) 
plot_pup_longum
ggsave("pups_longum.png", plot_pup_longum, dpi = 1200)
ggsave("pups_longum.svg", plot_pup_longum, device = "svg")

# no breve


######################dams
#####Phylum
###bars
ps_dam_phylum <- PSr_dam %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level 
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
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
ggsave("phylum_dams.png", plot6, dpi = 1200)
ggsave("phylum_dams.svg", plot6, device = "svg")

###Heatmap
ps_dam_phylum_z <- ps_dam_phylum %>%
  group_by(Phylum) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat6 <- ggplot(ps_dam_phylum_z, aes(x = diet, y = Phylum)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Phylum") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-1.7, 0.5)) +
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
ggsave("order_dams.png", plot7, dpi = 1200)
ggsave("order_dams.svg", plot7, device = "svg")

###heatmap
ps_dam_order_z <- ps_dam_order %>%
  group_by(Order) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat7 <- ggplot(ps_dam_order_z, aes(x = diet, y = Order)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Order") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-1.5, 0.5)) +
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
ggsave("family_dams.png", plot8, dpi = 1200)
ggsave("family_dams.svg", plot8, device = "svg")

###heatmap
ps_dam_family_z <- ps_dam_family %>%
  group_by(Family) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat8 <- ggplot(ps_dam_family_z, aes(x = diet, y = Family)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Family") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-1.5, 0.5)) +
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
ggsave("genus_dams.png", plot9, dpi = 1200)
ggsave("genus_dams.svg", plot9, device = "svg")

###heatmap
ps_dam_genus_z <- ps_dam_genus %>%
  group_by(Genus) %>%
  mutate(z_score = (Abundance - mean(Abundance)) / sd(Abundance)) %>%
  ungroup()

p.heat9 <- ggplot(ps_dam_genus_z, aes(x = diet, y = Genus)) + 
  geom_tile(aes(fill = z_score)) +
  labs(title = "Relative Abundance dams - Genus") +
  scale_fill_distiller("Relative\nAbundance\n(z-score)", palette = "RdYlBu", limits = c(-1.5, 0.5)) +
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
#getwd()
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/Bar_charts")
ggsave("species_dam.png", plot_dam_species, dpi = 1200)
ggsave("species_dam.svg", plot_dam_species, device = "svg")

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
ggsave("longum_dams.png", plot10, dpi = 1200)
ggsave("longum_dams.svg", plot10, device = "svg")

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
ggsave("longum_dams_zscore.png", plot10_z, dpi = 1200)
ggsave("longum_dams_zscore.svg", plot10_z, device = "svg")

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
ggsave("longum_dams_zscore_nolimits.png", plot10_z_nolimits, dpi = 1200)
ggsave("longum_dams_zscore_nolimits.svg", plot10_z_nolimits, device = "svg")

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
setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/Bar_charts")
ggsave("longum_dams_zscore_limits1025.png", plot10_z_limits1025, dpi = 1200)
ggsave("longum_dams_zscore_limits1025.svg", plot10_z_limits1025, device = "svg")

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
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/Bar_charts")
ggsave("longum_dams_zscore_limits310.png", plot10_z_limits310, dpi = 1200)
ggsave("longum_dams_zscore_limits310.svg", plot10_z_limits310, device = "svg")

#############stats for everything.Save as a table
#stat_pup_phylum <- ps_pup_phylum %>%
  #group_by(timepoint, sex) %>%
  #wilcox_test(Abundance ~ diet) %>%
  #add_xy_position(x = "diet")

#mutate(y.position = 100)
#stat_pup_phylum

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

setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance")
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

##subset w5M
#ps_pup_phylum_w5M <- ps_pup_phylum %>%
 # filter(timepoint == "w5" & sex == "M")

#stat_pup_phylum_w5M <- compare_means(data = ps_pup_phylum_w5M, 
 #                                    Abundance ~ diet) 
                                     #group.by = c("sex", "timepoint"))

##subset w10M
#ps_pup_phylum_w10M <- ps_pup_phylum %>%
 # filter(timepoint == "w10" & sex == "M")

#stat_pup_phylum_w10M <- compare_means(data = ps_pup_phylum_w10M, 
 #                                    Abundance ~ diet) 

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

setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/Abundance_tables_bysample")
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

############Generate abundance tables for each taxon by sample.
########Pups
phylum_abundance_table_sample_pup <- ps_pup_phylum %>%
  group_by(Phylum, Sample, diet, treatment, sex, timepoint, cohort, dam_cage, Abundance) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

order_abundance_table_sample_pup <- ps_pup_order %>%
  group_by(Order, Sample, diet, treatment, sex, timepoint, cohort, dam_cage, Abundance) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(sex, timepoint, diet, (Sample), Abundance)

family_abundance_table_sample_pup <- ps_pup_family %>%
  group_by(Family, Sample, diet, treatment, sex, timepoint, cohort, dam_cage, Abundance) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

genus_abundance_table_sample_pup <- ps_pup_genus %>%
  group_by(Genus, Sample, diet, treatment, sex, timepoint, cohort, dam_cage, Abundance) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

longum_abundance_table_sample_pup <- filtered_data %>%
  group_by(Genus_species, Sample, diet, treatment, sex, timepoint, cohort, dam_cage, Abundance) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(sex, timepoint, diet, desc(Abundance))

#setwd("C:/Users/Dara/OneDrive - University College Cork/ADATA onedrive/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance")
write_xlsx(phylum_abundance_table_sample_pup, "phylum_relativeabundance_sample_pup.xlsx")
write_xlsx(order_abundance_table_sample_pup, "order_relativeabundance_sample_pup.xlsx")
write_xlsx(family_abundance_table_sample_pup, "family_relativeabundance_sample_pup.xlsx")
write_xlsx(genus_abundance_table_sample_pup, "genus_relativeabundance_sample_pup.xlsx")
write_xlsx(longum_abundance_table_sample_pup, "longum_relativeabundance_sample_pup.xlsx")

########Dams
phylum_abundance_table_sample_dam <- ps_dam_phylum %>%
  group_by(Phylum, diet, Sample) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

order_abundance_table_sample_dam <- ps_dam_order %>%
  group_by(Order, diet, Sample) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

family_abundance_table_sample_dam <- ps_dam_family %>%
  group_by(Family, diet, Sample) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

genus_abundance_table_sample_dam <- ps_dam_genus %>%
  group_by(Genus, diet, Sample) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

longum_abundance_table_sample_dam <- filtered_datad %>%
  group_by(Genus_species, diet, Sample) %>% 
  summarize(Abundance = sum(Abundance, na.rm = TRUE)) %>% 
  arrange(diet, desc(Abundance))

#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance")
write_xlsx(phylum_abundance_table_sample_dam, "phylum_relativeabundance_sample_dam.xlsx")
write_xlsx(order_abundance_table_sample_dam, "order_relativeabundance_sample_dam.xlsx")
write_xlsx(family_abundance_table_sample_dam, "family_relativeabundance_sample_dam.xlsx")
write_xlsx(genus_abundance_table_sample_dam, "genus_relativeabundance_sample_dam.xlsx")
write_xlsx(longum_abundance_table_sample_dam, "longum_relativeabundance_sample_dam.xlsx")
