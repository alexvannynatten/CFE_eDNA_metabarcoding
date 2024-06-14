#####################################################################
# Last updated 2024-03-21 - Alexander Van Nynatten
# Compares the CFE mitofish amplicon sequences to the mitofish database
#####################################################################

#####################################################################
## Makes sequence database for new CFE sequences generated in this study
#####################################################################

# Inputs the new CFE sequences 
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ../00_raw_data/reference_sequences/saiab-cfe.fas \
--output-path temp/saiab-cfe.fas.qza

# Trims down to the mifish-U primer amplicon region
qiime feature-classifier extract-reads \
  --i-sequences temp/saiab-cfe.fas.qza \
  --p-f-primer GTCGGTAAAACTCGTGCCAGC \
  --p-r-primer CATAGTGGGGTATCTAATCCCAGTTTG \
  --p-min-length 100\
  --p-max-length 400\
  --o-reads temp/saiab-cfe.fas.trim.qza

#####################################################################
## Compares the sequence similarity of CFE MiFish in silico amplicons to mitofish database
#####################################################################

# Compares the reference amplicons to the mitofish database with vsearch
qiime feature-classifier classify-consensus-vsearch \
  --i-query temp/saiab-cfe.fas.trim.qza \
  --i-reference-reads only_mitofish_db.fas.trim.qza \
  --i-reference-taxonomy only_mitofish_db.tax.qza   \
  --p-maxaccepts 1000 \
  --p-top-hits-only \
  --o-classification temp/classified_rep_seqs_98.qza \
  --o-search-results temp/classified_blast6.qza

# outputs the results of the best match for each species in mitofish database
qiime tools export \
  --input-path temp/classified_blast6.qza \
  --output-path blast_results

#####################################################################