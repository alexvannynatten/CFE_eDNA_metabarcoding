#####################################################################
# Last updated 2024-03-21 - Alexander Van Nynatten
# Compares abundance calculated in mock community metabarcoding data to pooled DNA concentrations
#####################################################################

#####################################################################
## Loads libraries and files
#####################################################################

library(tidyverse)

ASV_table <- read_csv('04_asv_analyses/asv_sequence_table_corrected_onlyfish.csv')

#####################################################################
## Reformats the ASV table and subsamples mock community samples
#####################################################################

# Reformats ASV table to tidy dataframe
eDNA_df_tidy <- ASV_table %>%
  pivot_longer(!c(`OTU ID`, Fasta, o97, s97), 
               names_to = "sample.id", values_to = "Reads") %>%
  filter(sample.id %in% c('Mock_S1 Mock_S1', 'Mock_S2 Mock_S2')) %>%
  filter(Reads > 0) %>%
  separate_wider_delim(sample.id, " ", names = c("Location_code", "Extraction_code")) %>%
  filter(!s97 == 'unassigned') %>%
  rename(Species = s97) %>%
  group_by(Location_code, Species) %>%
  summarize(Reads = sum(Reads)) %>%
  ungroup() %>%
  group_by(Location_code) %>%
  mutate(Read_sums = sum(Reads)) %>%
  ungroup() %>%
  mutate(eDNA_frac = Reads / Read_sums) %>%
  ungroup() %>%

  # Pooled concentrations for each mock community
  mutate(Pool_frac = ifelse(Location_code == "Mock_S1", 1/6, 1/2)) %>%
  mutate(Difference = abs(eDNA_frac - Pool_frac))

#####################################################################
## Calculates how far off eDNA metabarcoding data is from expected values
#####################################################################

# Some summary stats
eDNA_df_tidy %>%
  ungroup() %>%
  summarize(Max_diff = max(Difference), 
            Min_diff = min(Difference), 
            Mean_diff = mean(Difference), 
            Median_diff = median(Difference))

#####################################################################
## Directly compares pooled DNA concentrations to DNA metabarcoding results
#####################################################################

# Formats dataframe for plotting
Fraction_plot <- eDNA_df_tidy %>%
  select(Species, Location_code, eDNA_frac, Pool_frac) %>%
  gather('Type', 'Fraction', -Location_code, -Species)

# Plots Figure S1
ggplot(data = Fraction_plot) + 
  geom_bar(
    aes(x = Fraction, y = Type, fill = Species), 
    colour = 'black', stat='identity', position = "fill"
  ) + 
  facet_grid(Location_code ~ .) + 
  theme_test() + 
  scale_fill_manual(values = c('#A75EC0', '#D53118', '#5DAF4C', '#287BC3', '#FFA243', '#D7D74B'))

ggsave('04_asv_analyses/FigS1_mock_community_plot.pdf')

#####################################################################