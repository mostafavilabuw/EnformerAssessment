This directory contins scripts that were used to make and analyze predictions for differential gene expression between 839 individuals
Please also refer to Fig1.sh, Fig2.sh, and FigS.sh to reproduce this analysis.

To generate the data used in these analyses, we use the Enformer model from 'https://tfhub.dev/deepmind/enformer/1', as described in 'https://github.com/deepmind/deepmind-research/tree/master/enformer'. 

# The following scripts use the Enformer model to generate attributions for differential gene expression prediction using the variants found in these 839 individuals. 
- These anlayses use variant information for each gene from ../variant_info_100k/{gene_id}.csv. These variant information files are produced from VCF file in 'EnformerAssessment/process_genomic_data/save_snp_info.txt'
-  These anlayses use reference sequences .../ref_seqs/{gene_id+}.npy' for each gene. These reference sequences are produced in 'EnformerAssessment/process_genomic_data/save_ref_seqs.py'

## per_SNP_ISM.py 
- save ISM results for Enformer on a given gene set by inserting the main (most common) variant for each SNP position in a gene 
- gene_file argument gives the genes to save this analysis for 

## gradient_attributions.py
- save (gradient at reference sequence x reference) and (gradient at reference sequence x main variant) for each main variant at each SNP position for a given gene
- gene_file argument gives the genes to save this analysis for 

## from_ref_drivers_ISM.py
- get ISM for each nucleotide within a window (here, 20bp) of each driver SNP (drivers found in select_drivers.py, as described below)
- start with the reference sequence, then insert each possible variant (same as from_main_var_drivers_ISM.py but without main variant inserted)
- edit 'drivers' in file to select which drivers to save this anlaysis for 

## from_main_var_drivers_ISM.py
- get ISM for each nucleotide within a window (here, 20bp) of each driver SNP (drivers found in select_drivers.py, as described below)
- start with the most common variant inserted, then insert each possible variant (same as from_ref_drivers_ISM.py but with main variant inserted)
- edit 'drivers' in file to select which drivers to save this anlaysis for 

## TSS_win_ISM.py 
- save Enformer output for each nucleotide inserted at every position within a window (here, 1000bp) of the TSS
- gene_file argument gives the genes to save this analysis for 

## predict_with_augmented_data.py
- test effects of data augmentation (+- 3bp shifts and reverse complement) on Enformer output
- gene_file argument gives the genes to save this analysis for 
- this analysis used individual genomes found in ../extractSequence/results/sequence100K/' (not available on github because of the size)



# The following scripts analyze the above attributions: 

## attribution_plot.py
- plot the attributions of all SNVs for a gene along the genomic location
`$python attribution_plot.py ENSG00000013573_ism_attributions.txt ISM --colors ../variant_info_100k/ENSG00000013573_frequency.txt --tss geneTSS.txt ENSG00000013573`

##  cluster_grad_attributions.py
-Clusters gradient attributions on the standard deviation within windows of 128bp
- Needs to be executed in directory with full length gradient and input sequence one-hot encodings (not available on github because of the size)
`$python cluster_grad_attributions.py <genesetfile> Prediction_correlationsCageAdultBrain.txt`

## compute_attribution.py
- compute attributions from prediction of variant and reference sequence
- execute in ism_res/ that contains <geneid>_ref_pred.npy
`$python compute_attribution.py <geneid>`

## compute_correlation.py
- compute correlation between observed and predicted expression across 839 individuals
`$python compute_correlation.py Observed_gene_expressionENSG.tsv Enformer_predictions_CageAdultBrain.txt > Prediction_correlationsCageAdultBrain.txt`
`$python compute_correlation.py Observed_gene_expressionENSG.tsv ../Enformer_predictions_CageAdultBrain.txt --get_allstats > Prediction_correlationsCageAdultBrain_Allstats.txt`

## compute_eqtl.py
- Compute eqtl values and correlation values for every SNVs genotype
`$python -W ignore compute_eqtl.py <genelistfile>`

## compute_population_frequency.py
- Compute the population frequency for every SNV
- execute in variant_info_100k/ with <geneid>snp_info.npz 
`$python compute_population_frequency.py <geneid>`

## compute_tstatpvalue.py 
- compute gene specific p-values for the predictions
`$python compute_tstatpvalue.py ../Prediction_correlationsCageAdultBrain.txt GeneSpecific_CorrelationtoObsRandomNull.txt`

## count_drivers.py
- Count the number of drivers per gene
- Execute in directory ism_res with files _ism_attributions_driversfw.txt
`$python count_drivers.py _ism_attributions_driversfw.txt > Counts_ism_attributions_driversfw.txt`

## driver_distance_to_tss.py
- Compute driver distance to TSS
`$python driver_distance_to_tss.py _ism_attributions_driversfw.txt geneTSS.txt > DistancetoTSS_ism_attributions_driversfw.txt`

## eqtl_attribution_plot.py
- Eqtl versus ISM plot for GSTM3
`$python eqtl_attribution_plot.py ENSG00000134202_ism_attributions.txt ../eqtl/ENSG00000134202_eqtl.txt ISM eQTL --colors ../variant_info_100k/ENSG00000134202_frequency.txt --drivers ENSG00000134202_ism_attributions_driversfw.txt --dpi 450 --fmt '.svg' --minfreq 0.01`

## eqtl_types.py
- Determine driver types
`$python eqtl_types.py <geneid>_ism_attributions.txt /eqtl/<geneid>_eqtl.txt <geneid>_ism_attributions_driversfw.txt`

## extract_ism_stats_around_drivers.py
- Analysis of motifs and attribution signals within and around drivers
- Execute in ism_drivers/, needs location of 'tss_attribution_stats.txt' and also access to /from_ref/ and /from_main_var/
`$python extract_ism_stats_around_drivers.py ../ism_res/ALLgenes_ism_attributions_driversfwmain_types.txt --checkmax --geneset <genesetfile>`

## find_common_motifs.py
- cluster one-hot encoded k-mers based on sequence similarity
`$python find_common_motifs.py one-hot_encoding.npz linkage(complete,single,average,etc.) cut-off(max distance for sequence to be included in cluster,f.e. 0.25 if 1/4 of bases can be distinct) min-length(minimum number of aligned bases between two sequences in a cluster) max-gap(maximum number of bases that don't match) --geneset <genesetfile>`

## generate_null.py
- Compute gene specific random null that is dependent on SNV structure
- Execute in variant_info_100k/ with <geneid>snp_info.npz
- Correct link to Observed gene expression file
`$python generate_null.py`

## plot_attribution_alongsequence.py
- Plot gradient for a gene along the entire sequence and also include drivers and other variants
`$python plot_attribution_alongsequence.py ENSG00000134202 0 ../ism_res/ENSG00000134202_ism_attributions.txt ../variant_info_100k/ENSG00000134202_frequency.txt ../ism_res/ENSG00000134202_ism_attributions_driversfw.txt ../ism_res/ENSG00000134202_ism_attributions_driversfw_types.txt --savefig`

## plot_common_motifs.py
- Check and plot enrichment of drivertypes for sequence clusters
`$python plot_common_motifs.py ALLgenes_ism_attributions_driversfw_types.txt ALLgenes_ism_attributions_driversfw_refseq_winsize13_clust_ms5-1_complete0.25 --savefig ALLgenes_ism_attributions_driversfw_refseq_winsize13_clust_ms5-1_complete0.2`

##  plot_distribution_enformer_correlations.py
- Plot violin plot for a selected set of genes versus the rest
`$python plot_distribution_enformer_correlations.py Prediction_correlationsCageAdultBrain.txt susie_SNP_gene_CortexENSG.txt`

## plot_drivercounts.py
- plot a histogram with the number of driver per gene
- Either plot for all genes or seperately for postively and negatively correlated genes
`$python plot_drivercounts.py Counts_ism_attributions_driversfw.txt ../Prediction_correlationsCageAdultBrain.txt ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --genelist ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list --split_sets`

## plot_driverdistance.py
- Plot histogram over the location of drivers with main drivers in different colors
`$python plot_driverdistance.py DistancetoTSS_ism_attributions_driversbw.txt ../Prediction_correlationsCageAdultBrain.txt ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --genelist ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list --plot_main_in_all`

##  plot_drivertype.py
- Plot enrichment of drivertypes for positive and negative correlated genes
`$python plot_drivertype.py ALLgenes_ism_attributions_driversfw_types.txt Prediction_correlationsCageAdultBrain.txt ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --genelist Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list -weighted ALLgenes_ism_attributions_driversfw.txt -3`

## plot_driver_motifstats.py
- Plot the ISM stats that were extracted with extract_ism_stats_around_drivers.py around the driver loci
`$python plot_driver_motifstats.py ../ism_res/ALLgenes_ism_attributions_driversfwmain_types.txt ALLgenes_ism_attributions_driversfwmain_ism_significance_stats.txt 2,15 --absolute --combine_columns max --scale 100 --lim 100 --cumulative -1 --nbins 22 --savefig ISMstats256_mainmaxzscore_in_var_or_ref.jpg`

## plot_individual_gene_scatter.py
- scatter plot between predicted and observed gene expression between individuals
`$python plot_individual_gene_scatter.py Observed_gene_expression.txt Enformer_predictions.txt DDX11 --figsize 4 3`

## plot_refandvar_attribution.py
- Generate ISM figures around a given loci
`$python ../plot_refandvar_attribution.py from_main_var/<geneid>_attribution.npz --norm_attributions ../tss_attribution_stats.txt <geneid> 2 --setylimpwm -1,1`
- Include also ISM of other SNVs in that area
- Mark drivers
- Include conservatin subplot
- exclude the per base effect heatmap
`$python plot_refandvar_attribution.py ENSG00000134202_109741038_attribution.npz --squaresize 0.4 --include_snvs 109741038 ../ism_res/ENSG00000134202_ism_attributions.txt ../eqtl/ENSG00000134202_corr.txt --markdrivers ../ism_res/ENSG00000134202_ism_attributions_driversfw.txt --dpi 350 --include_conservation ../PhyloP100/ENSG00000134202_1000tss_in_hg38.phyloP100way.txt --excludeheatmap --enlargepwm 1.8`

## plot_snp_clustering.py
- Plot distribution of number of variants per person for a gene
- Plot a sorted heatmap with SNVs as rows, individual genomes as columns and the genotype indicated by the color in the heatmap
`$python plot_snp_clustering.py ENSG00000013573snp_info.npz --minfreq 0.1 --minsnps 10 --combine_snps 0.9 --combine_individuals 0.9 --savefig ENSG00000013573_individual_enformer_snp`

## replace_genename.py 
- Replace gene names with ENSG-IDs in first column of a .txt file
`$python replace_genename.py Observed_gene_expression.txt`

## scatter_correlations.py
- Scatter plot comparison between PrediXcan and Enformer correlations to observed individual expression
`$python scatter_correlations.py --colors Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats.txt --colorscut 1.301 --vlim -1.5,1`

## scatter_mean_prediction.py 
- Plot mean expression of all individuals against predicted expression with Enformer from reference sequence 	
- change location of input file before using it
`$python scatter_mean_prediction.py`
* define test set
`$python scatter_mean_prediction.py --testset enformer_test_set_genes.npy`

## scatter_pvalue_vs_correlationprediction.py
- Plot p-value versus correlation of predicted expression with enformer from individual sequences
`$python scatter_pvalue_vs_correlationprediction.py Prediction_correlationsCageAdultBrain_Allstats.txt Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats.txt --markersize 4 --printset 0`

## select_drivers.py
- Determine drivers with forward method
- Correct link to file with predicted gene expression data
`$python select_drivers.py <geneid>_ism_attributions.txt ../variant_info_100k/<geneid>snp_info.txt <geneid> --forward`




