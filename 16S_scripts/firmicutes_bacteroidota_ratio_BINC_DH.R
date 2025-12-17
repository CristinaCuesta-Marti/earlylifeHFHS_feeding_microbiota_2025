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
library(ggsignif)
#install.packages("ggsignif")

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

##so removing blanks
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

#########################pups
####Firmicutes:Bacteroidetes ratio
PSr_pup_filtered <- PSr_pup %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt() %>%
  filter(Phylum %in% c("Firmicutes", "Bacteroidota")) %>%
  arrange(Phylum)

##Calculate the ratio pups
abundance_aggregated_pups <- PSr_pup_filtered %>%
  group_by(Sample, Phylum, diet, timepoint, sex) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = "drop")

abundance_wide_pups <- abundance_aggregated_pups %>%
  pivot_wider(names_from = Phylum, values_from = Total_Abundance, values_fill = 0) %>%
  mutate(ratio_pups = Firmicutes / Bacteroidota)

#write_xlsx(abundance_wide_pups, "pups_FBratio_table.xlsx")

###ratio fb plot
diet_colours <- c("Control Vehicle" = "#a0a0a4", "HFHS Vehicle" = "#f94040", "HFHS FOS+GOS" = "#addead", "HFHS B. longum APC1472" = "#e1c180")

plot_fb_pups_box <- ggplot(abundance_wide_pups, aes(x = diet, y = ratio_pups, fill = diet)) +
  geom_boxplot(outlier.shape = 16, position = position_dodge(width = 0.8), width = 0.7) +
  labs(x="", y= expression(italic("Firmicutes")*"/"*italic("Bacteroidota")*" ratio"), 
       title = expression(italic("Firmicutes")*"/"*italic("Bacteroidota")*" ratio in pups")) +
  facet_grid(sex ~ timepoint) +
  scale_fill_manual(values = diet_colours) +
  #geom_pwc(label = "p.adj.format", label.size = 3,
   #        p.adjust.method = "BH", vjust = 1.5, hide.ns = T, tip.length = 0.03) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        #legend.text = element_blank(),
        #legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot_fb_pups_box
#setwd("D:/BINC PROTECT study/Data_Analysis/Merged_before_assigning_taxonomy/maxEE24/Phyloseq_maxee24_items/Relative Abundance/Bacteroidota firmicutes/ratio")
ggsave("FB_ratio_pups_box_nopvalue.png", plot_fb_pups_box, dpi = 1200)
ggsave("FB_ratio_pups_box_nopvalue.svg", plot_fb_pups_box, device = "svg")

##ratio plot with p-values
plot_fb_pups_boxp <- ggplot(abundance_wide_pups, aes(x = diet, y = ratio_pups, fill = diet)) +
  geom_boxplot(outlier.shape = 16, position = position_dodge(width = 0.8), width = 0.7) +
  labs(x="", y= expression(italic("Firmicutes")*"/"*italic("Bacteroidota")*" ratio"), 
       title = expression(italic("Firmicutes")*"/"*italic("Bacteroidota")*" ratio in pups")) +
  facet_grid(sex ~ timepoint) +
  scale_fill_manual(values = diet_colours) +
  geom_pwc(label = "p.adj.format", label.size = 3,
          p.adjust.method = "BH", vjust = 1.5, hide.ns = T, tip.length = 0.03) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        #legend.text = element_blank(),
        #legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot_fb_pups_boxp

##export stats
p_values_pups <- abundance_wide_pups %>%
  group_by(sex, timepoint) %>%
  wilcox_test(ratio_pups ~ diet, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>% 
  mutate(
    group_comparison = paste(group1, "vs", group2), 
    formatted_p_value = format(p.adj, scientific = TRUE, digits = 3) 
  )
p_values_pups
write_xlsx(p_values_pups, "pups_FBratio_stats.xlsx")

######################Dams
PSr_dam_filtered <- PSr_dam %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt() %>%
  filter(Phylum %in% c("Firmicutes", "Bacteroidota")) %>%
  arrange(Phylum)

##Calculate the ratio. Method works
abundance_aggregated_dams <- PSr_dam_filtered %>%
  group_by(Sample, Phylum, diet) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = "drop")

abundance_wide_dams <- abundance_aggregated_dams %>%
  pivot_wider(names_from = Phylum, values_from = Total_Abundance, values_fill = 0) %>%
  mutate(ratio_dams = Firmicutes / Bacteroidota)
class(abundance_wide_dams)

#write_xlsx(abundance_wide_dams, "dams_FBratio_table.xlsx")

####ratio FB boxplot
plot_fb_dams_box <- ggplot(abundance_wide_dams, aes(x = diet, y = ratio_dams, fill = diet)) +
  geom_boxplot(outlier.shape = 16, position = position_dodge(width = 0.8), width = 0.7) +
  labs(x="", y= expression(italic("Firmicutes")*"/"*italic("Bacteroidota")*" ratio"), 
       title = expression(italic("Firmicutes")*"/"*italic("Bacteroidota")*" ratio in dams")) +
  #facet_grid(sex ~ timepoint) +
  scale_fill_manual(values = diet_colours) +
  #geom_pwc(label = "p.adj.format", label.size = 3,
   #        p.adjust.method = "BH", vjust = 1.8, hide.ns = T, tip.length = 0.03) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        #legend.text = element_blank(),
        #legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot_fb_dams_box
ggsave("FB_ratio_dams_box_nopvalue.png", plot_fb_pups_box, dpi = 1200)
ggsave("FB_ratio_pups_box_nopvalue.svg", plot_fb_pups_box, device = "svg")

###ratio FB with p values
plot_fb_dams_boxp <- ggplot(abundance_wide_dams, aes(x = diet, y = ratio_dams, fill = diet)) +
  geom_boxplot(outlier.shape = 16, position = position_dodge(width = 0.8), width = 0.7) +
  labs(x="", y= expression(italic("Firmicutes")*"/"*italic("Bacteroidota")*" ratio"), 
       title = expression(italic("Firmicutes")*"/"*italic("Bacteroidota")*" ratio in dams")) +
  #facet_grid(sex ~ timepoint) +
  scale_fill_manual(values = diet_colours) +
  geom_pwc(label = "p.adj.format", label.size = 3,
          p.adjust.method = "BH", vjust = 1.8, hide.ns = T, tip.length = 0.03) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(), 
        #legend.text = element_blank(),
        #legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot_fb_dams_boxp

##export stats
p_values_dams <- abundance_wide_dams %>%
  #group_by(sex, timepoint) %>%
  wilcox_test(ratio_dams ~ diet, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>% 
  mutate(
    group_comparison = paste(group1, "vs", group2), 
    formatted_p_value = format(p.adj, scientific = TRUE, digits = 3) 
  )
p_values_dams
write_xlsx(p_values_dams, "dams_FBratio_stats.xlsx")
