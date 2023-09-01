This directory contains intermediate data files that were generated with the scripts in enformer_analysis/ and can be used to reproduce intermediate results.

# ALLgenes_ism_attributions_driversfw_refseq_winsize13.npz
- Contains one-hot-encodings of 13bp around the driving SNVs
- The base in the center is the base that is present in the reference sequence

# ALLgenes_ism_attributions_driversfw_varseq_winsize13.npz
- Contains one-hot-encodings of 13bp around the driving SNVs
- The base in the center is the base is the variant instead of the reference

# ALLgenes_ism_attributions_driversfwmain_ism_significance_stats.txt
- Statistics from comparing the ISM values around the main driver of each gene with the maximum ISM that we observed within 2000bp around the TSS
- Almost as SupplementaryTableS3.txt

# ALLgenes_ism_attributions_driversfwmain_refseq_winsize13.npz
- Contains one-hot-encodings of 13bp around the main driving SNVs
- The base in the center is the base is the variant instead of the reference

# ALLgenes_ism_attributions_driversfwmain_varseq_winsize13.npz
- Contains one-hot-encodings of 13bp around the main driving SNVs
- The base in the center is the base is the variant instead of the reference

# Enformer_predictions.txt.gz
- Contains the predicted expression values for all expressed genes for 839 individual genotypes from the combination of all human output tracks with the fine-tuned elastic net model.

# MeanGeXPredFineTuned.txt
- Contains the predicted expression values for all expressed genes from their reference sequene from the combination of all human output tracks with the fine-tuned elastic net model.
- Contains the mean observed expression value

# Observed_gene_expression.txt.gz
- Contains the measured expression values for all expressed genes for 839 individual genotypes.

# PrediXcanCorrelationWithExpressionENSG.tsv
- Contains the Pearson correlation coefficient between the observed gene expression values of 839 individuals and the predicted values from PrediXcan.

# Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list
- Selected gene set to analyze with ISM
- Enformer correlation are significantly better than random models and absolute correlation is above 0.2

# Prediction_correlationsCageAdultBrain_Allstats.txt
- Contains Pearson correlation coefficients between the observed gene expression values of 839 individuals and the predicted values with Enformer's human CAGE,brain,adult track

# SupplementaryTable1.tsv
- Contains combined statistics from different files for all expressed genes

# SupplementaryTable2.tsv
- Contains statistics for all detected driver SNVs

# SupplementaryTable3.txt
- Contains statistics about the ISM values around all main drivers 

# enformer_test_set_genes.npy
- Gene set that was not used for training Enformer
- See Method section in paper for details

# gene-ids-and-positions.tsv
- Gene ids, gene names, location in hg38, location in hg19, and strand information 

# geneTSS.txt
- Transcription start sites

# ism_res.tar
- Contains directory ism_res/ which contains CAGE,brain,adult track prediction from reference sequence, and reference sequence with individually inserted main variants

# maindrivernpz.tar
- Contains .npz files with ISM values for 41 bp around the main drivers. 

# snp_positions.tar
- Contains the loci of all the SNVs that were within the window of Enformer inputs. 
- Loci match the variant predictions in ism_res.tar

# tss1000bpnpz.tar
- Contains .npz files with ISM values for +-1000 bp around the TSS

# tss_attribution_stats.txt
- Contains statistics from all files in tss1000bpnpz.tar, such as the max absolute ISM, and the standard deviation of ISMs within +-1000 bp around the TSS

# variant_info_100k.tar
- Personal gentype data cannot be shared with third party unless approved by the RADC. All requests to the investigator from the third party must be directed to the
RADC in order to assure knowledge of and compliance with the rules and regulations.
- Contains <geneid>snp_info.npz files
  - Files contain a matrix representing the genotype of all individuals, i.e rows are SNVs within input window and columns represent indivduals. 
  - Individuals with two copies of the major allele (genotype 0), those with one copy of the major allele (genotype 1) and those with two copies of the minor allele (genotype 2) 

