# Replace gene names with ENSG IDs
python3 replace_genename.py Observed_gene_expression.txt
python3 replace_genename.py Enformer_predictions.txt

# compute correlation between observed and predicted expression across 839 individuals
#python3 compute_correlation.py Observed_gene_expressionENSG.tsv ../Enformer_predictions_CageAdultBrain.txt > Prediction_correlationsCageAdultBrain.txt
#python3 compute_correlation.py Observed_gene_expressionENSG.tsv ../Enformer_predictions_CageAdultBrain.txt --get_allstats > Prediction_correlationsCageAdultBrain_Allstats.txt

# compute gene specific random null that is dependent on SNV structure
python3 generate_null.py

# compute gene specific p-values for the predictions
python3 compute_tstatpvalue.py ../Prediction_correlationsCageAdultBrain.txt GeneSpecific_CorrelationtoObsRandomNull.txt

# 2A: Plot p-value versus correlation of predicted expression with enformer from individual sequences
python3 scatter_pvalue_vs_correlationprediction.py Prediction_correlationsCageAdultBrain_Allstats.txt Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats.txt --markersize 4 --printset 0

# 2B: C2orf74 scatter plot
python3 plot_individual_gene_scatter.py Full_analysis/Observed_gene_expressionENSG.tsv Enformer_predictions_CageAdultBrain.txt ENSG00000237651
python3 plot_individual_gene_scatter.py Full_analysis/Observed_gene_expressionENSG.tsv Enformer_predictions_CageAdultBrain.txt ENSG00000134202 --setylim 0.057,0.071
python3 plot_individual_gene_scatter.py Full_analysis/Observed_gene_expressionENSG.tsv Enformer_predictions_CageAdultBrain.txt ENSG00000128944
python3 plot_individual_gene_scatter.py Full_analysis/Observed_gene_expressionENSG.tsv Enformer_predictions_CageAdultBrain.txt ENSG00000133433
python3 plot_individual_gene_scatter.py Full_analysis/Observed_gene_expressionENSG.tsv Enformer_predictions_CageAdultBrain.txt ENSG00000226752
python3 plot_individual_gene_scatter.py Full_analysis/Observed_gene_expressionENSG.tsv Enformer_predictions_CageAdultBrain.txt ENSG00000120675

# 2C: Scatter plot comparison between PrediXcan and Enformer correlations to observed individual expression
python3 scatter_correlations.py --colors Full_analysis/Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats.txt --colorscut 1.301 --vlim -1.5,1

# compute eqtl
cd eqtl/
python3 -W ignore compute_eqtl.py ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list
cd ..

cd ism_res/
# determine drivers with forward method
genes=$(ls *_ism_attributions.txt)
for g in genes
do
python3 ../select_drivers.py $g ../variant_info_100k/${g%_ism_attributions.txt}snp_info.txt ${g%_ism_attributions.txt} --forward
done 

# 2D: Eqtl versus ISM plot for GSTM3
python3 ../eqtl_attribution_plot.py ENSG00000134202_ism_attributions.txt ../eqtl/ENSG00000134202_eqtl.txt ISM eQTL --colors ../variant_info_100k/ENSG00000134202_frequency.txt --drivers ENSG00000134202_ism_attributions_driversfw.txt --dpi 450 --fmt '.svg' --minfreq 0.01

# determine driver types
genes=$(ls *_ism_attributions.txt)
for g in genes
do
python3 ../eqtl_types.py $g ../eqtl/${g%_ism_attributions.txt}_eqtl.txt ${g%_ism_attributions.txt}_ism_attributions_driversfw.txt
done

# combine driver types
python3 ../combine_eqtl_types.py _ism_attributions_driversfw_types.txt > ALLgenes_ism_attributions_driversfw_types.txt

# 2E: Enrichment of drivertypes for positive and negative correlated genes
python3 ../plot_drivertype.py ALLgenes_ism_attributions_driversfw_types.txt ../Prediction_correlationsCageAdultBrain.txt ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --genelist ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list -weighted ALLgenes_ism_attributions_driversfw.txt -3

# Compute driver distance to TSS
python3 ../driver_distance_to_tss.py _ism_attributions_driversfw.txt ../geneTSS.txt > DistancetoTSS_ism_attributions_driversfw.txt

# 2F: Location of drivers with main drivers in different colors
python3 ../plot_driverdistance.py DistancetoTSS_ism_attributions_driversbw.txt ../Prediction_correlationsCageAdultBrain.txt ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --genelist ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list --plot_main_in_all

# 2G: Gradient plot for GSTM3 with drivers and other variants
python3 ../plot_attribution_alongsequence.py ENSG00000134202 0 ../ism_res/ENSG00000134202_ism_attributions.txt ../variant_info_100k/ENSG00000134202_frequency.txt ../ism_res/ENSG00000134202_ism_attributions_driversfw.txt ../ism_res/ENSG00000134202_ism_attributions_driversfw_types.txt --savefig

# Generate ISM figure around the TSS
python3 ../plot_refandvar_attribution.py ENSG00000134202_109741038_attribution.npz --squaresize 0.12 --include_snvs 109741038 ../ism_res/ENSG00000134202_ism_attributions.txt ../eqtl/ENSG00000134202_corr.txt --markdrivers ../ism_res/ENSG00000134202_ism_attributions_driversfw.txt --dpi 350 --include_conservation ../PhyloP100/ENSG00000134202_1000tss_in_hg38.phyloP100way.txt --excludeheatmap --enlargepwm 1.8
python3 ../plot_refandvar_attribution.py $g --figsize 75,1 --include_snvs $loc ../ism_res/${gene}_ism_attributions.txt ../eqtl/${gene}_corr.txt --markdrivers ../ism_res/${gene}_ism_attributions_driversfw.txt --dpi 250 --include_conservation ../PhyloP100/${gene}_1000tss_in_hg38.phyloP100way.txt --excludeheatmap --enlargepwm 1.2

