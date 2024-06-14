#####################################################################
# Last updated 2024-03-21 - Alexander Van Nynatten
# Uses CUTADAPT and DADA2 (Qiime2) to process metabarcoding data
#####################################################################

#####################################################################
## Establishes the environment for the analysis
#####################################################################

#####################################################################
## Generates statistics for the quality of raw sequencing data

# Imports paired end reads from fastq files
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../../00_raw_data/sample_metadata/sample_manifest.txt \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path ../../03_qiime2_pipeline/temp/MiSeq_raw_reads.qza

cd ../../03_qiime2_pipeline/temp

# Visualizes the statistics regarding the raw reads
qiime demux summarize \
  --i-data MiSeq_raw_reads.qza \
  --o-visualization MiSeq_raw_reads.qzv

#####################################################################
## Trims primers and heterogeneity spacers with CUTADAPT in Qiime2
#####################################################################

# Removes primers and heterogeniety spacers
qiime cutadapt trim-paired \
--i-demultiplexed-sequences MiSeq_raw_reads.qza \
  --p-front-f GTCGGTAAAACTCGTGCCAGC \
  --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
  --p-error-rate 0.2 \
  --p-discard-untrimmed \
  --o-trimmed-sequences MiSeq_trimmed_reads.qza

# Visualizes the statistics regarding the raw reads
qiime demux summarize \
  --i-data MiSeq_trimmed_reads.qza \
  --o-visualization MiSeq_trimmed_reads.qzv

#####################################################################
## Denoising and dereplicating reads using DADA2
#####################################################################

# Denoises reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs MiSeq_trimmed_reads.qza \
  --p-trunc-len-f 110 \
  --p-trunc-len-r 110 \
  --o-table dada2_table.qza \
  --o-representative-sequences dada2_rep_set.qza \
  --o-denoising-stats dada2_stats.qza

# Outputs list of the frequency of reads
qiime feature-table summarize \
  --i-table dada2_table.qza \
  --m-sample-metadata-file ../../00_raw_data/sample_metadata/sample_metadata.txt \
  --o-visualization dada2_table.qzv

# Outputs list of representative sequences
qiime feature-table tabulate-seqs \
  --i-data dada2_rep_set.qza \
  --o-visualization rep-seqs.qzv

# Outputs table of the denoised sequences
qiime metadata tabulate \
  --m-input-file dada2_stats.qza \
  --o-visualization dada2_stats.qzv

#####################################################################