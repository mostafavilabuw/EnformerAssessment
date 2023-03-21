# Files and scripts to generate the supplementary figures

# S1. Prediction of mean expression for the test set
python3 scatter_mean_prediction.py --testset enformer_test_set_genes.npy

#S2. Distribution of number of variants per person for DDX11
python3 ../plot_snp_clustering.py ENSG00000013573snp_info.npz --select_individuals ../Enformer_predictions_individuals.txt --minfreq 0.1 --minsnps 10 --combine_snps 0.9 --combine_individuals 0.9 --savefig ENSG00000013573_individual_enformer_snp

#S4 Scatter plot between correlations to observed for fine-tuned and CAGE-track
python3 scatterplot.py Prediction_correlationsCageAdultBrain_Allstats.txt ism_res/Refpred.txt "CAGE,brain,adult,MeanSumlog10+-2indv" "CAGE,brain,adult,log10sum+-1ref" --columns -2 -1 --density --alpha 0.5 --label --filternan --linewidth 0.1 --log10y

#NOT included: Mean of Predicted CAGE tracks and sum of Observed, both are logged
python3 scatterplot.py Prediction_correlationsCageAdultBrain_Allstats.txt Prediction_correlationsCageAdultBrain_Allstats.txt "MeanObs" "CAGE,brain,adult,MeanSumlog10+-2indv" --columns -4 -2 --density --alpha 0.5 --label --filternan --linewidth 0.1
python3 scatterplot.py Prediction_correlationsCageAdultBrain_Allstats.txt Prediction_correlationsCageAdultBrain_Allstats.txt "StdObs" "CAGE,brain,adult,StdSumlog10+-2indv" --columns -3 -1 --density --alpha 0.5 --label --filternan --linewidth 0.1
python3 scatterplot.py Prediction_correlationsCageAdultBrain_Allstats.txt Prediction_correlationsCageAdultBrain_Allstats.txt "MeanObs" "StdObs" --columns -4 -3 --density --alpha 0.5 --label --filternan --linewidth 0.1
python3 scatterplot.py Prediction_correlationsCageAdultBrain_Allstats.txt Prediction_correlationsCageAdultBrain_Allstats.txt "MeanCAGE,adult,brain" "StdCAGE,adult,brain" --columns -2 -1 --density --alpha 0.5 --label --filternan --linewidth 0.1

#S6 Distribution of Susie genes
python3 plot_distribution_enformer_correlations.py Full_analysis/Prediction_correlationsCageAdultBrain.txt susie_SNP_gene_CortexENSG.txt

# S7. Correlation between ISM and gradient attributions between variants of all investigated genes
python3 correlation_ism_grad_attributions.py

# compute the sum of snp attributions
cd ism_res/
genes=$(ls ENSG*_ism_attributions.txt)
for g in $genes
do
python3 ../sum_meanattributions.py $g ../variant_info_100k/${g%_ism_attributions.txt}snp_info.npz
done
# combine the sum of predictions
python3 ../combine_predictions.py ENSG _ism_attributions_sumpersonal_mp.txt 
# compute correlation to Enformers predictions
python3 ../compute_correlation.py ALL_genes_ism_attributions_sumpersonal_mp.txt ../Enformer_predictions_CageAdultBrain.txt > ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt
cd ..

# compute sum of gradient attributions
cd ref_attribs/
genes=$(ls *_grad_attributions.txt)
for g in $genes
do
time python3 ../sum_meanattributions.py $g ../variant_info_100k/${g%_grad_attributions.txt}snp_info.npz
done
# combine the sum of predictions
python3 ../combine_predictions.py ENSG _grad_attributions_sumpersonal_mp.txt
# compute correlation to Enformers predictions
python3 ../compute_correlation.py ALL_genes_grad_attributions_sumpersonal_mp.txt ../Enformer_predictions_CageAdultBrain.txt > ALL_genes_grad_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt
cd ..

# S8. Correlation between ISM and gradient sums for every individual with predicted values of enformer for every individual
python3 plot_correlations_predictions_tosumattributions.py Prediction_correlationsCageAdultBrain.txt ism_res/ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt ref_attribs/ALL_genes_grad_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --savefig

# Not shown anymore: Intermediate plots for forward identification of drivers
# GSTM3
python3 ../select_drivers.py ENSG00000134202_ism_attributions.txt ../variant_info_100k/ENSG00000134202snp_info.txt ENSG00000134202 --forward --plot_test 40
# DDX11
python3 ../select_drivers.py ENSG00000013573_ism_attributions.txt ../variant_info_100k/ENSG00000013573snp_info.txt ENSG0000013573 --forward --plot_test 40

# combine drivers
python3 ../combine_drivers.py _ism_attributions_driversfw.txt > ALLgenes_ism_attributions_driversfw.txt

# S9. Plot number of SNP drivers
python3 ../count_drivers.py _ism_attributions_driversfw.txt > Counts_ism_attributions_driversfw.txt
python3 ../count_drivers.py _ism_attributions_driversbw.txt > Counts_ism_attributions_driversbw.txt

python3 ../plot_drivercounts.py Counts_ism_attributions_driversfw.txt ../Prediction_correlationsCageAdultBrain.txt ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --genelist ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list --split_sets

# S10 sequence patterns at the driver snps

#Make list with main drivers
python3 ../select_main_drivers.py ALLgenes_ism_attributions_driversfw.txt > ALLgenes_ism_attributions_driversfwmain.txt
# make type file with main drivers
python3 ../select_driverset.py ALLgenes_ism_attributions_driversfwmain.txt ALLgenes_ism_attributions_driversfw_types.txt ALLgenes_ism_attributions_driversfwmain_types.txt

# extract sequence windows of size 13 around drivers
python3 ../extract_sequence_around_drivers.py ALLgenes_ism_attributions_driversfw.list 6 ../geneTSS.txt
python3 ../extract_sequence_around_drivers.py ALLgenes_ism_attributions_driversfwmain.list 6 ../geneTSS.txt
# cluster 13-mers for ref seq
python3 ../find_common_motifs.py ALLgenes_ism_attributions_driversfw_refseq_winsize13.npz complete 0.25 5 1 --geneset ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list
python3 ../find_common_motifs.py ALLgenes_ism_attributions_driversfwmain_refseq_winsize13.npz complete 0.25 5 1 --geneset ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list
# check enrichment of specific types of clusters
python3 ../plot_common_motifs.py ALLgenes_ism_attributions_driversfw_types.txt ALLgenes_ism_attributions_driversfw_refseq_winsize13_clust_ms5-1_complete0.25 --savefig ALLgenes_ism_attributions_driversfw_refseq_winsize13_clust_ms5-1_complete0.2 
# check enrichment only for main drivers
python3 ../plot_common_motifs.py ALLgenes_ism_attributions_driversfwmain_types.txt ALLgenes_ism_attributions_driversfwmain_refseq_winsize13_clust_ms5-1_complete0.25 --savefig ALLgenes_ism_attributions_driversfwmain_refseq_winsize13_clust_ms5-1_complete0.25

# cluster 13mers for varseq
python3 ../find_common_motifs.py ALLgenes_ism_attributions_driversfw_varseq_winsize13.npz complete 0.25 5 1 --geneset ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list
python3 ../find_common_motifs.py ALLgenes_ism_attributions_driversfwmain_varseq_winsize13.npz complete 0.25 5 1 --geneset ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list
# check enrichment of specific types of clusters
python3 ../plot_common_motifs.py ALLgenes_ism_attributions_driversfw_types.txt ALLgenes_ism_attributions_driversfw_varseq_winsize13_clust_ms5-1_complete0.25 --savefig ALLgenes_ism_attributions_driversfw_varseq_winsize13_clust_ms5-1_complete0.25
# check enrichment only for main drivers
python3 ../plot_common_motifs.py ALLgenes_ism_attributions_driversfwmain_types.txt ALLgenes_ism_attributions_driversfwmain_varseq_winsize13_clust_ms5-1_complete0.25 --savefig ALLgenes_ism_attributions_driversfwmain_varseq_winsize13_clust_ms5-1_complete0.25

# S11. Full size figures of driver distance from tss
# see Fig2.sh

# S12. Clusters of gradient attribution locations
python3 ../cluster_grad_attributions.py ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list ../Prediction_correlationsCageAdultBrain.txt

#Not included. Causal SNP enrichment in positive genes
#python3 susie_SNP_gene_analysis.py susie_SNP_gene_brain_union_b.txt --genelist Full_analysis/Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.txt
## potentially 
#python3 susie_SNP_gene_analysis.py susie_SNP_gene_CortexENSG.txt --genelist Full_analysis/Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.txt


# S13. Analysis of motifs and attribution signals within and around drivers
# Analyze signal at drivers from ism
python3 extract_ism_stats_around_drivers.py ../ism_res/ALLgenes_ism_attributions_driversfwmain_types.txt --checkmax --geneset Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list

python3 ../plot_driver_motifstats.py ../ism_res/ALLgenes_ism_attributions_driversfwmain_types.txt ALLgenes_ism_attributions_driversfwmain_ism_significance_stats.txt 2,15 --absolute --combine_columns max --scale 100 --lim 100 --cumulative -1 --nbins 22 --savefig ISMstats256_mainmaxzscore_in_var_or_ref.jpg

python3 ../plot_driver_motifstats.py ../ism_res/ALLgenes_ism_attributions_driversfwmain_types.txt ALLgenes_ism_attributions_driversfwmain_ism_significance_stats.txt 4,17 --combine_columns max --lim 10 --cumulative -1 --nbins 11 --savefig ISMstats256_mainmotifssize_in_var_or_ref_10percofmax.jpg --xlim 0,10

python3 ../plot_driver_motifstats.py ../ism_res/ALLgenes_ism_attributions_driversfwmain_types.txt ALLgenes_ism_attributions_driversfwmain_ism_significance_stats.txt 5,18 --combine_columns max --lim 10 --cumulative -1 --nbins 11 --savefig ISMstats256_mainmotifssize_in_var_or_ref_20percofmax.jpg --xlim 0,10

python3 ../plot_driver_motifstats.py ../ism_res/ALLgenes_ism_attributions_driversfwmain_types.txt ALLgenes_ism_attributions_driversfwmain_ism_significance_stats.txt 6,19 --combine_columns max --lim 10 --cumulative -1 --nbins 11 --savefig ISMstats256_mainmotifssize_in_var_or_ref_50percofmax.jpg --xlim 0,10

genes='ENSG00000001460_24417389 ENSG00000013573_31073901 ENSG00000134202_109741163'
for g in $genes
do
gene=${g::15}
python3 ../plot_refandvar_attribution.py from_main_var/${g}_attribution.npz --norm_attributions ../TSS_ISM/tss_attribution_stats.txt ${gene} 2 --setylimpwm -1,1
done
genes='ENSG00000001460_24416934_1000bp_attribution.npz ENSG00000013573_31073845_1000bp_attribution.npz ENSG00000134202_109741038_1000bp_attribution.npz'
for g in $genes
do
gene=${g::15}
loc=${g%_1000bp_attribution.npz}
loc=${loc:16}
python3 ../plot_refandvar_attribution.py $g --norm_attributions  --norm_attributions ../TSS_ISM/tss_attribution_stats.txt ${gene} 2 --figsize 100,3 --include_snvs $loc ../ism_res/${gene}_ism_attributions.txt ../variant_info_100k/${gene}_frequency.txt --markdrivers ../ism_res/${gene}_ism_attributions_driversfw.txt --dpi 150 --include_conservation ../PhyloP100/${gene}_1000tss_in_hg38.phyloP100way.txt --excludeheatmap
done


