Environmental DNA metabarcoding in the Cape Fold aquatic ecoregion: opportunities and challenges for eDNA uptake in an endemism hotspot 
=====================

**Objective of study**

To evaluate the performance of eDNA metabarcoding in highly endemic fish communities and compare the detection and relative abundance estimates obtained from this data to conventional underwater visual surveys.

**Overview or dataset**
eDNA metabarcoding data and visual survey data from eight sampling sites in the Cape Fold aquatic ecoregion (CFE). Reference 12S sequences for endemic fishes.

--------------------------

**00_raw_data/**: Directory containing all raw files for analyses
- **map_data/** 
	- Shape files available from: https://www.dws.gov.za/iwqs/gis_data/river/rivs500k.aspx
- **miseq_fastq/** 
	- Paired-end Illumina reads archived SRAXXXXXXX. In format:
		- *R1_001.fastq.gz*
		- *R2_001.fastq.gz*
- **reference_sequences/** 
	- *mitofish-all.fas*: full mitofish database, version 3.87 Update: 2023-03-08, available from: https://mitofish.aori.u.tokyo.ac.jp
	- *saiab-cfe.fas*: new 12S sequences generated in this study for fishes from the CFE.
- **sample_metadata/** 
	- *location_metadata.csv*: eDNA sampling / visual survey locations (lat/long) and environmental metadata.
	- *reference_metadata.csv*: location data for museum specimens used for 12S reference database (lat/long).
	- *sample_manifest.txt*: Qiime2 formatted manifest file for MiSeq paired-end reads.
	- *sample_metadata.txt*: Qiime2 formatted metadata for each MiSeq sample.
- **visual_surveys/** 
	- *visual_counts_output.csv*: Raw visual survey data (snorkel and underwater camera)
	- *visual_species_codes.csv*: Database used to correct simplified species names in *visual_counts_output.csv*

--------------------------

*01_1_CFE_map.R*: Plots the Cape Fold Eco-region boundaries, eDNA sampling locations, and the locations where reference tissue samples were collected for sequence database. 

**01_sample_map/**: Directory containing all output for figure 1 (Map of CFE and sampling locations).
- *Fig1a_Africa.pdf*: Map of Africa, CFE highlighted.
- *Fig1a_CFE.pdf*: Detailed map of CFE.
- *Fig1bc_DrainageMap.pdf*: Map of sampling locations in Bos and Blindekloof Rivers.

--------------------------

*02_1_qiime2_database.sh*: Formats the mitofish and cfe fasta files into Qiime2 sequence and taxonomy files for ASV classification using the ```feature-classifier extract-reads``` and ```tools import``` function in qiime2.
*02_2_qiime2_blast.sh*: Compares the CFE mitofish amplicon sequences to the mitofish-only reference database using the ```feature-classifier classify-consensus-vsearch``` function in qiime2.

**02_reference_database/**
	- **blast_results/**: results for Table S4.
	- *cfe_mitofish_db.fas.trim.qza*: Qiime2 sequence artifact.
	- *cfe_mitofish_db.tax.qza* Qiime2 taxonomy artifact.
	- *only_mitofish_db.fas.trim.qza* Qiime2 sequence artifact.
	- *only_mitofish_db.tax.qza* Qiime2 taxonomy artifact.
	- **temp_data/**: Directory where intermediate qiime2 artifacts are stored for qiime2_database pipeline.

--------------------------

*03_1_qiime2_denoising.sh*: Qiime2 pipeline to remove primer sequences with ```qiime cutadapt trim-paired``` and merge/denoise/dereplicate amplicons using ```qiime dada2 denoise-paired```.  
*03_2_qiime2_classification.sh*: Qiime2 pipeline used to classify sequences using ```qiime feature-classifier classify-consensus-vsearch```. 
*03_3_asv_tidying.R*: Corrects for contamination using read counts from NTCs and removes low frequency reads.

**03_qiime2_pipeline/**
	- *asv_sequence_table_o97.tsv*: The uncorrected ASV table classified using the mitofish-only reference databases.
	- *asv_sequence_table_s97.tsv*: The uncorrected ASV table classified using the mitofish reference databases supplemented with the new 12S sequences from CFE fishes.
	- *asv_sequence_table_corrected.csv*: An ASV table where sequencing reads are corrected for contamination using controls and low frequency reads are removed. Contains classification information for both datasets.
	- **temp_data/**

--------------------------

*04_1_classification_comparison.R*: Compares the number of read counts classified between reference databases by rarifying read counts across samples using ```phyloseq::rarefy_even_depth```.
*04_2_mock_communities.R*: Compares abundance calculated in mock community metabarcoding data to pooled DNA concentrations.
*04_3_insilico_test.R*: In silico PCR comparison of the efficiency of two 12S primer sets on 12S sequences from CFE fishes.

**04_asv_analyses/**

- *asv_sequence_table_corrected_onlyfish.csv*: ASV table where all non-fish sequences are removed. Used in all subsequent analyses.
- *Fig2_sp_classification.pdf*: Comparison of the proportion of eDNA metabarcoding sequencing reads.
- *FigS1_mock_community_plot.pdf*: Comparison of relative read abundance in mock communities to the DNA concentrations pooled.
- *in_silico_PCR_df.csv*: Table of PCR efficiencies calculated by *04_3_insilico_test.R*

--------------------------

*05_1_detection_analyses.R*: Compares species detections at sites from eDNA and visual methods using a contingency table and calculates the coefficient of correlation using ```psych::phi```.
*05_2_abundance_analyses.R*: Compares species detections at sites from eDNA and visual methods. Correlation between each method was tested using ```stats::lm```.

**05_statistical_analyses/**

- *Fig3_sp_detections.pdf*: Comparison of detections of taxa by eDNA metabarcoding and visual surveys.
- *Fig4a_sp_abundance.pdf*: Comparison of eDNA metabarcoding relative read frequency to abundance data from Underwater cameras.
- *Fig4b_sp_abundance.pdf*: Comparison of eDNA metabarcoding relative read frequency to abundance data from Snorkel surveys.
