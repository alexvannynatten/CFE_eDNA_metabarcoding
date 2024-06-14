#####################################################################
# Last updated 2023-03-21 - Alexander Van Nynatten
# Corrects for contamination using NTCs and removes low frequency reads
#####################################################################

#####################################################################
# Loads files and libraries 
#####################################################################

library(tidyverse)

# Loads all qiime2 output files
data_files <- list.files('./03_qiime2_pipeline', 
                         pattern=".tsv", full.names = TRUE)

metadata_df <- read.table('./00_raw_data/sample_metadata/sample_metadata.txt', 
                          sep = "\t", header = TRUE)

#####################################################################
## Makes tidy dataframe of the qiime2 output asv tables
#####################################################################

# Transforms the qiime2 taxonomy string to separate columns
combined_df <- data.frame()
for(i in 1:length(data_files)) {
  
  ASV_table <- read_delim(data_files[[i]], delim = '\t')
  
  ASV_seq_df <- ASV_table %>%
    pivot_longer(!c(`OTU ID`, Taxon, Consensus, Fasta), 
                 names_to = "sample.id", values_to = "Count") %>%
    filter(Count > 0) %>%
    mutate(Taxon = gsub('D_[a-z]__', '', Taxon)) %>%
    separate_wider_delim(Taxon, ";", names = c("order",
                                               "family", 
                                               "genus", 
                                               "species"), 
                         too_few = "align_start") %>%
    mutate(Database = gsub('.*_', '', data_files[[i]])) %>%
    mutate(Database = gsub("\\.tsv$", "", Database))
  
  # Combines output for both databases
  combined_df <- rbind(combined_df, ASV_seq_df)
  
}

#####################################################################
# Adds metadata for each sample
#####################################################################

uncorrected_df <- combined_df %>%
  select(!c(order, family, genus, Consensus)) %>%
  replace_na(list(species = 'unassigned')) %>%
  pivot_wider(names_from = Database, values_from = species) %>%
  left_join(metadata_df) %>%
  ungroup()

#####################################################################
# Removes contaminant reads based on blank filters and NTCs
#####################################################################

# Gets the max count for each OTU in blank or NTC sample
Max_contam <- uncorrected_df %>%
  filter(Sample_type %in% c('Blank', 'NTC')) %>%
  group_by(`OTU ID`) %>%
  summarize(max_control = max(Count)) %>%
  ungroup()

corrected_df <- uncorrected_df %>%

  # Removes reads based on the maximum number observed for each OTU in the control samples
  left_join(Max_contam) %>%
  rowwise() %>%
  replace_na(list(max_control = 0)) %>%
  ungroup() %>%
  mutate(Count = Count - max_control) %>%
  filter(Count > 0) %>%

  # Removes short reads
  filter(nchar(Fasta) > 150) %>%

  # Pools reads across sample PCR replicates
  group_by(`OTU ID`, o97, s97, Location_code, Extraction_code, Fasta, Sample_type) %>%
  summarize(Count = sum(Count)) %>%
  ungroup() %>%
  group_by(Extraction_code) %>%
  mutate(Count_sum = sum(Count)) %>%
  ungroup() %>%
  mutate(Count_perc = Count / Count_sum) %>%

  # Removes low frequency ASVs
  filter(Count_perc > 0.013)

#####################################################################
# Converts corrected dataframes back into ASV table formats 
#####################################################################

# Generates the ASV table
asv_table_output <- corrected_df %>%
  ungroup() %>%
  arrange(Extraction_code) %>%
  mutate(Extraction_id = paste(Location_code, Extraction_code)) %>%
  select(`OTU ID`, Fasta, o97, s97, Count, Extraction_id) %>%
  pivot_wider(names_from = Extraction_id, values_from = Count) %>%
  replace(is.na(.), 0)

write.csv(asv_table_output, '03_qiime2_pipeline/asv_sequence_table_corrected.csv', row.names = FALSE, quote = FALSE)

#####################################################################
