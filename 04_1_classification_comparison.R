#####################################################################
# Last updated 2024-03-21 - Alexander Van Nynatten
# Compares read counts between databases by rarifying reads
#####################################################################

#####################################################################
# Loads files and libraries 
#####################################################################

set.seed(01)

library(tidyverse)
library(phyloseq)

ASV_table <- read_csv('04_asv_analyses/asv_sequence_table_corrected_onlyfish.csv') %>%
  select(!c('Mock_S1 Mock_S1', 'Mock_S2 Mock_S2'))

#####################################################################
# Rarifies read count across samples 
#####################################################################

# Makes matrix from ASV table
asvmat <- apply(as.matrix.noquote(ASV_table[ ,c(5:ncol(ASV_table))]),2,as.numeric)
rownames(asvmat) <- ASV_table$`OTU ID`

# Makes phyloseq asvmat formatted file
ASV = otu_table(asvmat, taxa_are_rows = TRUE)

# Rarifies reads
ASV.rarified <- rarefy_even_depth(ASV, sample.size = min(sample_sums(ASV)),
                                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Tidy dataframe from the rarified ASV matrix
ASV_rarefied_output <- data.frame(ASV.rarified) %>%
  mutate(`OTU ID` = row.names(.)) %>%
  pivot_longer(-`OTU ID`, names_to = 'Extraction_id', values_to = 'Count') %>%
  filter(Count > 0) %>%
  left_join(ASV_table[ ,c(1:4)]) %>%
  separate_wider_delim(Extraction_id, ".", names = c("Location_code", "Extraction_code")) %>%
  mutate(River = ifelse(grepl('BK', Location_code), 'BK', 'Bos'))

#####################################################################
# Plots the comparison of classifications between the two datasets
#####################################################################

Classification_db <- ASV_rarefied_output %>%
  select(River, Count, s97, o97) %>%
  pivot_longer(names_to = 'Database', cols = ends_with("97"), values_to = 'Species') %>%
  group_by(Species, River, Database) %>%
  summarize(Count = sum(Count))

# Plots Figure 2
ggplot(data = Classification_db) + 
  geom_bar(
    aes(x = River, y = Count, fill = Species), 
    colour = 'black', stat='identity', position = "fill"
  ) + 
  facet_grid(. ~ Database) + 
  theme_test() + 
  scale_fill_manual(values = c('#A75EC0', 
                               '#01A3E0', 
                               '#D53118', 
                               '#5DAF4C', 
                               '#F55433', 
                               '#16BCCE', 
                               '#287BC3', 
                               '#FFA243', 
                               '#D7D74B', 
                               '#D4D4D4'))

ggsave('04_asv_analyses/Fig2_sp_classification.pdf')

#####################################################################
# Calculates the proportion of reads classified with either database
#####################################################################

# Total reads generated for fishes
Total_rare_reads <- sum(Classification_db$Count)/2

# Fraction of reads unclassified in each reference database
Classification_db %>%
  group_by(Database, Species) %>%
  summarize(Count_total = sum(Count)) %>%
  filter(Species == 'unassigned') %>%
  mutate(1 - (Total_unclassified = Count_total / Total_rare_reads))


#####################################################################
# Number of sequences per sample and site
#####################################################################

# Some summary data for the sequencing analysis
ASV_table %>%
  pivot_longer(!c(`OTU ID`, Fasta, o97, s97), 
               names_to = "sample.id", values_to = "Reads") %>%
  filter(Reads > 0) %>%
  separate_wider_delim(sample.id, " ", names = c("Location_code", "Extraction_code")) %>%
  group_by(Location_code) %>%
  summarize(Read_sums = sum(Reads)) %>%
  ungroup() %>%
  summarize(avg_reads = mean(Read_sums))

#####################################################################
  