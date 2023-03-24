This directory contains the paramaters and Matlab script that was used to predict fine-tuned gene expression using the Enformer model.

The Enformer model used can be found in https://github.com/deepmind/deepmind-research/tree/master/enformer. 
We separately predict gene expression for the paternal and maternal sequence for each individual. 
The script enformerFinetuned.m is used to finetune the predicted gene expression by assigning a weight to each of Enformer's 5313 output tracks. 
