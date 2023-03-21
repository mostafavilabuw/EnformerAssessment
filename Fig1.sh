# 1C: Plot mean expression of all individuals against predicted expression with enformer from reference sequence
python3 scatter_mean_prediction.py

# 1D: DDX11 scatter plot 
python3 plot_individual_gene_scatter.py Observed_gene_expression.txt Enformer_predictions.txt DDX11 --figsize 4 3

# List of all genes with npy files to get attributions
gene=$(ls *_ref_pred.npy)
# Compute attributions all genes
for g in $gene
do
echo ${g%_ref_pred.npy}
python3 compute_attribution.py ${g%_ref_pred.npy}
done

gene=$(ls *snp_info.npz)
# Compute the population frequency for every SNV
for g in $gene
do
python3 ../compute_population_frequency.py ${g%snp_info.txt}
done

# 1E: DDX11 ISM attribution plot
python3 ../attribution_plot.py ENSG00000013573_ism_attributions.txt ISM --colors ../variant_info_100k/ENSG00000013573_frequency.txt --tss ../geneTSS.txt ENSG00000013573


