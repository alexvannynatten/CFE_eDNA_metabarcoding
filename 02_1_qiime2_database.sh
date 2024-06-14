#####################################################################
# Last updated 2024-03-21 - Alexander Van Nynatten
# Formats the fasta and taxonomic files for Qiime2 classification
#####################################################################

#####################################################################
## Generates qiime2 refseq and tax files for the mitofish database exclusively
#####################################################################

cat ../00_raw_data/reference_sequences/*all.fas >> temp/only_mitofish_db.fas

# Reformats the sequence, order, family, genus, species (seq, ord, fam, gen, spp) names for the taxonomy file
grep ">" temp/only_mitofish_db.fas | awk -F '[ ]' '{print $1}' | sed 's/>//' > temp/o_seq
grep ">" temp/only_mitofish_db.fas | awk -F '[ ]' '{print $2}' | sed 's/[^A-Za-z]*//g' | sed 's;^;\tD_o__;' > temp/o_ord
grep ">" temp/only_mitofish_db.fas | awk -F '[ ]' '{print $3}' | sed 's/[^A-Za-z]*//g' | sed 's/^\s*$/NA/' | sed 's/^/;D_f__/' > temp/o_fam
grep ">" temp/only_mitofish_db.fas | awk -F '[ ]' '{print $1}' | awk -F '[|]' '{print $3}' | awk -F '[_]' '{print $1}' | sed 's/^/;D_g__/' > temp/o_gen
grep ">" temp/only_mitofish_db.fas | awk -F '[ ]' '{print $1}' | awk -F '[|]' '{print $3}' | tr -s '_' ' ' | sed 's/^/;D_s__/' > temp/o_spp

# Concatenates taxonomy information for the taxonomy file
paste -d'\0' temp/o_seq temp/o_ord temp/o_fam temp/o_gen temp/o_spp > temp/only_mitofish_db.tax

# Reformats the sequence names for fasta file
cut -d ' ' -f 1 < temp/only_mitofish_db.fas > temp/only_mitofish_db.fa

# Imports taxonomy strings corresponding to sequence file
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path temp/only_mitofish_db.fa \
--output-path temp/only_mitofish_db.fas.qza

# Trims down to the mifish-U primer amplicon region
qiime feature-classifier extract-reads \
  --i-sequences temp/only_mitofish_db.fas.qza \
  --p-f-primer GTCGGTAAAACTCGTGCCAGC \
  --p-r-primer CATAGTGGGGTATCTAATCCCAGTTTG \
  --p-min-length 100\
  --p-max-length 400\
  --o-reads only_mitofish_db.fas.trim.qza

# Imports taxonomy strings corresponding to sequence file
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path temp/only_mitofish_db.tax \
  --output-path only_mitofish_db.tax.qza

#####################################################################
## Generates qiime2 refseq and tax files for the combined mitofish and new cfe 12S sequences database
#####################################################################

cat ../00_raw_data/reference_sequences/*.fas >> temp/cfe_mitofish_db.fas

# Reformats the sequence, order, family, genus, species (seq, ord, fam, gen, spp) names for the taxonomy file
grep ">" temp/cfe_mitofish_db.fas | awk -F '[ ]' '{print $1}' | sed 's/>//' > temp/s_seq
grep ">" temp/cfe_mitofish_db.fas | awk -F '[ ]' '{print $2}' | sed 's/[^A-Za-z]*//g' | sed 's;^;\tD_o__;' > temp/s_ord
grep ">" temp/cfe_mitofish_db.fas | awk -F '[ ]' '{print $3}' | sed 's/[^A-Za-z]*//g' | sed 's/^\s*$/NA/' | sed 's/^/;D_f__/' > temp/s_fam
grep ">" temp/cfe_mitofish_db.fas | awk -F '[ ]' '{print $1}' | awk -F '[|]' '{print $3}' | awk -F '[_]' '{print $1}' | sed 's/^/;D_g__/' > temp/s_gen
grep ">" temp/cfe_mitofish_db.fas | awk -F '[ ]' '{print $1}' | awk -F '[|]' '{print $3}' | tr -s '_' ' ' | sed 's/^/;D_s__/' > temp/s_spp

# Concatenates taxonomy information for the taxonomy file
paste -d'\0' temp/s_seq temp/s_ord temp/s_fam temp/s_gen temp/s_spp > temp/cfe_mitofish_db.tax

# Reformats the sequence names for fasta file
cut -d ' ' -f 1 < temp/cfe_mitofish_db.fas > temp/cfe_mitofish_db.fa

# Imports taxonomy strings corresponding to sequence file
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path temp/cfe_mitofish_db.fa \
--output-path temp/cfe_mitofish_db.fas.qza

# Trims down to the mifish-U primer amplicon region
qiime feature-classifier extract-reads \
  --i-sequences temp/cfe_mitofish_db.fas.qza \
  --p-f-primer GTCGGTAAAACTCGTGCCAGC \
  --p-r-primer CATAGTGGGGTATCTAATCCCAGTTTG \
  --p-min-length 100\
  --p-max-length 400\
  --o-reads cfe_mitofish_db.fas.trim.qza

# Imports taxonomy strings corresponding to sequence file
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path temp/cfe_mitofish_db.tax \
  --output-path cfe_mitofish_db.tax.qza

#####################################################################
