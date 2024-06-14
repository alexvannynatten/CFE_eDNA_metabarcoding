#####################################################################
# Last updated 2024-03-21 - Alexander Van Nynatten
# Compares species detections at sites from eDNA and visual methods
#####################################################################

#####################################################################
## Loads the data and libraries required
#####################################################################

library(tidyverse)
library(psych)

# Loads site and sequencing read data
samples_df <- read_csv('00_raw_data/sample_metadata/location_metadata.csv')
sequences_df <- read_csv('04_asv_analyses/asv_sequence_table_corrected_onlyfish.csv')
visual_df <- read_csv('00_raw_data/visual_surveys/visual_counts_output.csv')
sp_codes <- read.csv('00_raw_data/visual_surveys/visual_species_codes.csv')

#####################################################################
# Generates dataframe of the visual detections for plotting
#####################################################################

# Orders river by latitude to match river image
samples_df <- arrange(samples_df, lat) %>%
  mutate(Location_code = factor(Location_code, 
                                levels = unique(Location_code)))

# Makes list of detections by visual survey methods 
visual_df_tidy <- visual_df %>%
  select(-c(2,3)) %>%
  gather('Method', 'Count', -Location_code) %>%
  mutate(Sp_code = gsub("_.*","",Method),
         Method = sub("_[^_]+$", "",Method)) %>%
  mutate(Method = gsub(".*_","",Method)) %>%
  filter(!Method == 'density') %>%
  left_join(sp_codes) %>%
  filter(Count > 0) %>%
  select('Location_code', 'Species', 'Method') %>%
  unique() %>%
  mutate(n = 1) %>%
  mutate(Location_code = factor(Location_code, 
                                levels = levels(samples_df$Location_code))) %>%
  mutate(River = ifelse(grepl('BK', Location_code), 'BK', 'Bos')) %>%
  mutate(Pseudobarbus_sp = ifelse(River == 'BK', 'afer', 'swartzi')) %>%
  mutate(Species = ifelse(Species == 'Pseudobarbus', paste(Species, Pseudobarbus_sp), Species)) %>%
  select(-Pseudobarbus_sp, -River)

#####################################################################
# Generates dataframe of the eDNA based detections for plotting
#####################################################################

eDNA_df_tidy <- sequences_df %>%
  pivot_longer(!c(`OTU ID`, Fasta, o97, s97), 
               names_to = "sample.id", values_to = "Count") %>%
  filter(Count > 0) %>%
  separate_wider_delim(sample.id, " ", names = c("Location_code", "Extraction_code")) %>%
  filter(!Location_code %in% c('Mock_S1', 'Mock_S2')) %>%
  filter(!s97 == 'unassigned') %>%
  group_by(Location_code, s97) %>%
  summarize(Reads = sum(Count), n = n_distinct(Extraction_code)) %>%
  ungroup() %>%
  select(Location_code, s97, n) %>%
  mutate(Method = 'eDNA') %>%
  na.omit() %>%
  rename(Species = s97)
  
#####################################################################
# Plots the comparison of the two methods
#####################################################################

# Plots figure 3
ggplot() + 
  geom_point(data = visual_df_tidy,
             aes(x= Species, y = Location_code),
             size = 12, shape = 22, colour = 'black', fill = 'grey90') + 
  geom_point(data = eDNA_df_tidy,
             aes(x= Species, y = Location_code, fill = Species),
             shape = 21, size = 6, colour = 'black') + 
  geom_text(data = eDNA_df_tidy,
            aes(x= Species, y = Location_code, label = n)) + 
  theme_test() + 
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_fill_manual(values = c('#A75EC0', '#01A3E0', '#D53118', '#5DAF4C', '#F55433', '#16BCCE', '#287BC3', '#FFA243', '#D7D74B', '#D4D4D4'))

ggsave('05_statistical_analyses/Fig3_sp_detections.pdf')

#####################################################################
## Statistically compares the consistency of detections between both methods
#####################################################################

# Makes contingency plot
Contingency_table <- rbind(visual_df_tidy, eDNA_df_tidy) %>%
  mutate(River = ifelse(grepl('BK', Location_code), 'BK', 'Bos')) %>%
  mutate(Method = ifelse(Method %in% c('camera', 'snorkle'), 'visual', Method)) %>%
  distinct() %>%
  mutate(n = n^0) %>%
  group_by(River) %>%
  complete(Species, Method, nesting(Location_code), fill = list(n = 0)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = Method, 
    values_from = n) %>%
  select(-Location_code, -Species, -River) %>%
  table()

# Calculates the stats
phi(Contingency_table)

#####################################################################