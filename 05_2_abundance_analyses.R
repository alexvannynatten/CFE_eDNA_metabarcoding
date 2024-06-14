#####################################################################
# Last updated 2024-03-21 - Alexander Van Nynatten
# Compares relative abundance of species at sites using eDNA and visual methods
#####################################################################

#####################################################################
## Loads the data and libraries required for analyses
#####################################################################

library(tidyverse)
library(lmtest)

# Loads site and sequencing read data
sequences_df <- read_csv('04_asv_analyses/asv_sequence_table_corrected_onlyfish.csv')
visual_df <- read_csv('00_raw_data/visual_surveys/visual_counts_output.csv')
sp_codes <- read.csv('00_raw_data/visual_surveys/visual_species_codes.csv')

#####################################################################
## Generates dataframe of the visual detections for plotting
#####################################################################

# Makes list of detections by visual survey methods 
vis_df_tidy <- visual_df %>%
  select(-c(2,3)) %>%
  gather('Type', 'Count', -Location_code) %>%
  mutate(Sp_code = gsub("_.*","",Type),
         Type = gsub(".*_","",Type)) %>%
  left_join(sp_codes) %>%
  mutate(River = ifelse(grepl('BK', Location_code), 'BK', 'Bos')) %>%
  mutate(Pseudobarbus_sp = ifelse(River == 'BK', 'afer', 'swartzi')) %>%
  mutate(Species = ifelse(Species == 'Pseudobarbus', paste(Species, Pseudobarbus_sp), Species)) %>%
  select(-Pseudobarbus_sp, -River) %>%
  select('Location_code', 'Type', 'Count', 'Species') %>%
  filter(Type %in% c('maxN', 'average')) %>%
  group_by(Location_code, Type) %>%
  mutate(Sum = sum(Count)) %>%
  mutate(Vis_frac = Count / Sum)

#####################################################################
# Generates dataframe of the eDNA based detections for plotting
#####################################################################

eDNA_df_tidy <- sequences_df %>%
  pivot_longer(!c(`OTU ID`, Fasta, o97, s97), 
               names_to = "sample.id", values_to = "Reads") %>%
  filter(Reads > 0) %>%
  separate_wider_delim(sample.id, " ", names = c("Location_code", "Extraction_code")) %>%
  filter(!Location_code %in% c('Mock_S1', 'Mock_S2')) %>%
  filter(!s97 == 'unassigned') %>%
  rename(Species = s97) %>%
  mutate(River = ifelse(grepl('BK', Location_code), 'BK', 'Bos')) %>%
  group_by(Location_code, Extraction_code, Species, River) %>%
  summarize(Reads = sum(Reads)) %>%
  ungroup() %>%
  group_by(Location_code, Extraction_code) %>%
  mutate(Read_sums = sum(Reads)) %>%
  ungroup() %>%
  mutate(eDNA_frac = Reads / Read_sums) %>%
  ungroup()

#####################################################################
## Merges the visual abundance and eDNA-based abudnance values into a single dataframe
#####################################################################

# Makes new dataframe with both types of data
# change "maxN" to "average" for figure 4b
Comparison_df <- eDNA_df_tidy %>%
  left_join(vis_df_tidy %>%
              filter(Type == 'maxN')
  ) %>%

  # Removes and instances where only one species was detected
  filter(eDNA_frac > 0) %>%
  filter(Vis_frac > 0) %>%
  mutate(Big_pool = ifelse(Location_code %in% c('BK17', 'BK18'), 'Big', 'Small')) %>%
  filter(!Location_code %in% c('B14', 'B17')) %>%
  data.frame()

# Removes the two largest pools surveyed
Comparison_df_small <- Comparison_df %>%
  filter(!Big_pool == 'Big')

#####################################################################
## Calculates the variance between separate eDNA replicates
#####################################################################

group_variances <- tapply(Comparison_df$eDNA_frac, paste(Comparison_df$Location_code, Comparison_df$Species), var)
median(group_variances[!is.na(group_variances)])

#####################################################################
## Tests the linear relationship between the visual and eDNA abundace
#####################################################################

# Linear mode
lm_model <- lm(eDNA_frac ~ Vis_frac, data = Comparison_df)
summary(lm_model)

# Residuals
Comparison_df$Residuals <- lm_model$residuals

# Boxplot of residuals
ggplot(Comparison_df) + 
  geom_boxplot(aes(x = Big_pool, y = abs(Residuals))) + 
  theme_test()

# Wilcox test comparing residuals for large and small pools
wilcox.test(abs(Residuals) ~ Big_pool, data = Comparison_df)

#####################################################################
## Tests the linear relationship between the visual and eDNA abundace for only small pools
#####################################################################

lm_model <- lm(eDNA_frac ~ Vis_frac, data = Comparison_df_small)
summary(lm_model)

# Checks assumptions of linear models are met
shapiro.test(residuals(lm_model))
qqnorm(residuals(lm_model))
qqline(residuals(lm_model))
plot(lm_model, which = 1)
bptest(lm_model)

#####################################################################
## Plots the data
#####################################################################

# Plots Figure 4
ggplot(data = Comparison_df,
       aes(x = Vis_frac, y = eDNA_frac)) + 
  geom_smooth(aes(group = Big_pool),
              method = "lm", se = TRUE) +
  geom_point(aes(fill = as.character(Species), colour = as.character(Species),
                 shape = Big_pool),
             size = 3) + 
  theme_test() + 
  scale_shape_manual(values = c(3, 21)) +
  xlim(0,1) +
  ylim(0,1) + 
  scale_fill_manual(values = c('#D53118', '#5DAF4C', '#16BCCE', '#287BC3', '#FFA243', '#D7D74B')) +
  scale_colour_manual(values = c('#D53118', '#5DAF4C', '#16BCCE', '#287BC3', '#FFA243', '#D7D74B'))

ggsave('05_statistical_analyses/Fig4a_sp_abundance.pdf')

#####################################################################

