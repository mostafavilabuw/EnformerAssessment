%% Finetune Enformer output with GTEx mean expression model
% Input:    enformPat = #subjects x 5313 tracks enformer output matrix of a gene for paternal sequences
%           enformMat = #subjects x 5313 tracks enformer output matrix of a gene for maternal sequences
%           modelPath = Path where GTEx finetuned model is stored
% Output:   exprPred = #subjects x 1 predicted gene expression levels
% Note:     Assumes enformPat and enformMat = log2 of enformer outputs AFTER summing over output windows e.g. +/-2 bins from center
function exprPred = enformerFinetuned(enformPat,enformMat,modelPath)
%% Load mean expression model
fid = fopen(modelPath);
header = regexp(fgetl(fid),',','split');
temp = textscan(fid,'%s%f','Delimiter',',');
fclose(fid);
wt = temp{2};
beta = wt(1:end-1);
a0 = wt(end);

%% Enformer outputs for Paternal
exprPredPat = enformPat*beta+a0;

%% Enformer outputs for Maternal
exprPredMat = enformMat*beta+a0;

%% Combine paternal and maternal effect
exprPred = exprPredPat+exprPredMat;
