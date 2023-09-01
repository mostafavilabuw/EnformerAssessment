# How far are we from personalized gene expression prediction using sequence-to-expression deep neural networks?

This repository contains the data and the scripts that were used to assess Enformer's [[1]](#1) ability to predict differential gene expression across 839 individuals [[2]](#2).

To reproduce our analysis [[2]](#2), please follow the examplary commands in Fig1.sh, Fig2.sh, and FigS.sh. These bash scripts were not intended to run by themselves. Instead, these scripts were meant to guide interested people through our analysis and document our analysis. Please **run each step separately** and if necessary replace hard coded file locations in python scripts with the correct location in the Data directory. Most intermediate output files from these scripts can also be found in the Data/ directory, so that intermediate steps can be skipped. 

To reproduce the entire analysis, it is **essential to obtain access to personal genotype data from the the Rush Memory and Aging Project (ROSMAP)**. The personal genotype data that was used in this study cannot be shared with a third party unless approved by the RADC. All requests to the investigator from the third party must be directed to the RADC in order to assure knowledge of and compliance with the rules and regulations. Genotype, RNA-seq, and DNAm data for the Religious Orders Study and Rush Memory and Aging Project (ROSMAP) samples are available from the Synapse AMP-AD Data Portal (Accession Code: syn2580853) as well as RADC Research Resource Sharing Hub at https://www.radc.rush.edu/.

For RNA-seq pre-processing, please refer to the Supplementary Methods in [[2]](#2). In brief, we applied TMM normalization (using edgeR calcNormFactors) to the raw counts to estimate the effective library size of each individual. We then applied voom/limma to regress out confounds and convert the counts into log2(CPM). Technical covariates included: batch, study (ROS or MAP), RNA integrity number, postmortem interval, Library size, log PF number of aligned reads, PCT_CODING_BASES, PCT_INTERGENIC_BASES, PCT_PF_READS_ALIGNED, PCT_RIBOSOMAL_BASES, PCT_UTR_BASES, PERCENT_DUPLICATION, MEDIAN_3PRIME_BIAS, MEDIAN_5PRIME_TO_3PRIME_BIAS,  MEDIAN_CV_COVERAGE. Biological covariates, including 1) age, 2) sex, and 3) top 10 expression principal components. Both biological and technical covariates were regressed out from log raw read counts. Only genes with mean log2(CPM) > 2 were included. Mean expression values were retained for downstream analysis.

The variant call files for whole genome sequencing (WGS) data from the ROSMAP in variant call format (VCF) were obtained from the Synapse repository (syn117074200). The coordinates of variant calls (GRCh37) were converted to GRCh38 coordinates using the Picard LiftoverVcf tool (http://broadinstitute.github.io/picard). The Eagle software2 version 2.4.1 was used to phase the genotypes with the default setting.

## Abstract

_Deep learning (DL) methods accurately predict various functional properties from genomic DNA, including gene expression, promising to serve as an important tool in interpreting the full spectrum of genetic variations in personal genomes. However, systematic out-of-sample benchmarking is needed to assess the gap in their utility as personalized DNA interpreters. Using paired Whole Genome Sequencing and gene expression data we evaluate DL sequence-to-expression models, identifying their critical failure to make correct predictions on a substantial number of genomic loci, highlighting the limits of the current model training paradigm._

## References
<a id="1">[1]</a> 
Avsec, Ž., Agarwal, V., Visentin, D., Ledsam, J.R., Grabska-Barwinska, A., Taylor, K.R., Assael, Y., Jumper, J., Kohli, P. and Kelley, D.R., 2021. Effective gene expression prediction from sequence by integrating long-range interactions. [Nature methods](https://www.nature.com/articles/s41592-021-01252-x), 18(10), pp.1196-1203.

<a id="2">[2]</a>
Sasse, A.\*, Ng, B.\*, Spiro, A.\*, Tasaki, S., Bennett, D., Gaiteri, C., De Jager, P.L., Chikina, M., and Mostafavi, S., 2023. How far are we from personalized gene expression prediction using sequence-to-expression deep neural networks? [bioRxiv](https://doi.org/10.1101/2023.03.16.532969), https://doi.org/10.1101/2023.03.16.532969, \* These authors contributed equally
