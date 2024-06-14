#####################################################################
# Last updated 2024-03-21 - Alexander Van Nynatten
# Classifies sequences using VSERACH (qiime2) and exports data
#####################################################################

#####################################################################
## Classifies sequences using only mitofish database
#####################################################################

# Classifies representative sequences to species using VSEARCH
qiime feature-classifier classify-consensus-vsearch \
  --i-query dada2_rep_set.qza \
  --i-reference-reads ../../02_reference_database/only_mitofish_db.fas.trim.qza \
  --i-reference-taxonomy ../../02_reference_database/only_mitofish_db.tax.qza \
  --p-perc-identity 0.97 \
  --p-top-hits-only \
  --p-maxaccepts 100 \
  --o-classification classified_rep_seqs_o97.qza \
  --o-search-results classified_blast6_o97.qza

# Outputs a list of the classifications given to each sequence
qiime metadata tabulate \
  --m-input-file classified_rep_seqs_o97.qza \
  --o-visualization classified_rep_seqs_o97.qzv

# Visualization of the number of classified reads in each sample
qiime taxa barplot \
  --i-table dada2_table.qza \
  --i-taxonomy ./classified_rep_seqs_o97.qza \
  --m-metadata-file ../../00_raw_data/sample_metadata/sample_metadata.txt \
  --o-visualization ./taxa_barplot_o97.qzv

# Exports classified reads to csv files based on taxonomic level
qiime tools export \
  --input-path classified_rep_seqs_o97.qzv \
  --output-path class_table_o97

#####################################################################
## Classifies sequences using CFE supplemented mitofish database
#####################################################################

# Classifies representative sequences to species using VSEARCH
qiime feature-classifier classify-consensus-vsearch \
  --i-query dada2_rep_set.qza \
  --i-reference-reads ../../02_reference_database/cfe_mitofish_db.fas.trim.qza \
  --i-reference-taxonomy ../../02_reference_database/cfe_mitofish_db.tax.qza \
  --p-perc-identity 0.97 \
  --p-top-hits-only \
  --p-maxaccepts 100 \
  --o-classification classified_rep_seqs_s97.qza \
  --o-search-results classified_blast6_o97.qza

# Outputs a list of the classifications given to each sequence
qiime metadata tabulate \
  --m-input-file classified_rep_seqs_s97.qza \
  --o-visualization classified_rep_seqs_s97.qzv

# Visualization of the number of classified reads in each sample
qiime taxa barplot \
  --i-table dada2_table.qza \
  --i-taxonomy ./classified_rep_seqs_s97.qza \
  --m-metadata-file ../../00_raw_data/sample_metadata/sample_metadata.txt \
  --o-visualization ./taxa_barplot_s97.qzv

# Exports classified reads to csv files based on taxonomic level
qiime tools export \
  --input-path classified_rep_seqs_s97.qzv \
  --output-path class_table_s97

#####################################################################

# Exports reads counts for each sequence in each sample
qiime tools export \
--input-path dada2_table.qza \
--output-path read_table

biom convert \
-i ./read_table/feature-table.biom \
-o ./read_table/asv_table.tsv --to-tsv

# Exports classified reads to csv files based on taxonomic level
qiime tools export \
  --input-path rep-seqs.qzv \
  --output-path seq_table

#####################################################################
## Generates final ASV tables combining sequence classifications, read tables, and amplicon fastas
#####################################################################

tail -n +3 class_table_o97/metadata.tsv | sort > class_table_o97/qiime2_class_seqs_98.tsv.sort
echo -e "$(echo -n 'OTU ID\tTaxon\tConsensus\n'; cat class_table_o97/qiime2_class_seqs_98.tsv.sort)" > class_table_o97/qiime2_class_seqs_98.tsv.sort

tail -n +3 class_table_s97/metadata.tsv | sort > class_table_s97/qiime2_class_seqs_98.tsv.sort
echo -e "$(echo -n 'OTU ID\tTaxon\tConsensus\n'; cat class_table_s97/qiime2_class_seqs_98.tsv.sort)" > class_table_s97/qiime2_class_seqs_98.tsv.sort

# need -e if bash (if zsh remove it)
tail -n +3 read_table/asv_table.tsv | sort > read_table/asv_table.tsv.sort
echo -e "$(echo -n "$(sed -n '2p' read_table/asv_table.tsv | cut -c2-)\n"; cat read_table/asv_table.tsv.sort)" > read_table/asv_table.tsv.sort

awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1"\t"$2}' seq_table/sequences.fasta | sort > seq_table/qiime2_sequences.fasta.sort
echo -e "$(echo -n 'OTU ID\tFasta\n'; cat seq_table/qiime2_sequences.fasta.sort)" > seq_table/qiime2_sequences.fasta.sort

# Combines files into one sequence table containing read counts, classifcations, and fasta sequences
join -t $'\t' -j 1 read_table/asv_table.tsv.sort class_table_o97/qiime2_class_seqs_98.tsv.sort | join -t $'\t' -j 1 - seq_table/qiime2_sequences.fasta.sort | tr -d "\015" > ../asv_sequence_table_o97.tsv
join -t $'\t' -j 1 read_table/asv_table.tsv.sort class_table_s97/qiime2_class_seqs_98.tsv.sort | join -t $'\t' -j 1 - seq_table/qiime2_sequences.fasta.sort | tr -d "\015" > ../asv_sequence_table_s97.tsv

#####################################################################