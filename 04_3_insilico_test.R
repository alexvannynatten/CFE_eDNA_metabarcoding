#####################################################################
# Last updated 2024-03-21 - Alexander Van Nynatten
# In silico PCR comparison of two different 12S primer sets for fishes
#####################################################################

# Libraries
library(DECIPHER)
library(tidyverse)

# CFE fishes 12S sequences
Seq_db <- readDNAStringSet('00_raw_data/reference_sequences/saiab-cfe.fas')

# Simplifies names for analysis
Simple_names <- data.frame(Names = names(Seq_db)) %>%
  mutate(Names = gsub('.*\\|', '', Names)) %>%
  mutate(Names = gsub('_.*', '', Names))
  
names(Seq_db) <- Simple_names$Names

#####################################################################
## Loops through primer sequences and reference sequence testing primer efficiency
#####################################################################

# Primer sequences
primer_pairs <- list(
mifish = c('GTCGGTAAAACTCGTGCCAGC', 'CATAGTGGGGTATCTAATCCCAGTTTG'),
Teleo = c('ACACCGCCCGTCACTCT', 'CTTCCGGTACACTTACCATG')
)

# In silico PCR of all sequences in the database with primers above
is_PCR_df <- data.frame()
for(i in 1:length(primer_pairs)){
  for(j in 1:length(Seq_db)){

is_PCR <- AmplifyDNA(primer_pairs[[i]],
        Seq_db[j],
        maxProductSize = 1000,
        annealingTemp = 52,
        P = 4e-7,
        ions = 0.2,
        includePrimers=TRUE,
        )

if(length(is_PCR) < 1){
next
}

# Formats output
is_PCR_df <- rbind(is_PCR_df, 
  data.frame(Species = sub("_[^_]+$", "", names(Seq_db[j])), 
             Sequence = as.character(is_PCR[[1]]), 
             Efficiency = max(as.numeric(gsub('%.*', '', names(is_PCR)))),
             Primers = names(primer_pairs[i])
             )
  )
print(j)
}
}

#####################################################################
## Summarizes the primer efficiencies for each species
#####################################################################

# Calculates the median efficiency for each primer, 0 if not amplified
is_PCR_df_sum <- is_PCR_df %>%
  expand(Species, Primers) %>%
  rowwise() %>%
  mutate(Efficiency = median(as.numeric(is_PCR_df$Efficiency)[
    paste(is_PCR_df$Species, is_PCR_df$Primers) %in% 
    paste(Species, Primers)])) %>%
  replace_na(list(Efficiency = 0))

write.csv(is_PCR_df_sum, '04_asv_analyses/in_silico_PCR_df.csv', row.names=FALSE)

#####################################################################