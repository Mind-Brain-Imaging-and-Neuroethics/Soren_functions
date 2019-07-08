function p = Permtest_ISC(data1,data2,nperm,type)
% Permtest_ISC uses a permutation test to test for difference in means 
% between two intersubject correlation matrices
%
% Input arguments: 
%
% data1 and data2 are the data from which the ISC matrices are defined. 
%    They should be in the form observations x subjects
% nperm is the number of permutations to use
%
% type should be either 'pearson' or 'spearman'
%
% this method assumes the same number of observations in each group
%
% Output arguments:
%
% p is the p-value from the permutation test


if nargin < 4
    type = 'pearson';
end

corrmat1 = corr(data1,'Type',type);
corrmat2 = corr(data2,'Type',type);
corrvals1 = corrmat1(find(belowDiag(ones(size(corrmat1)))));
corrvals2 = corrmat2(find(belowDiag(ones(size(corrmat2)))));
[~,~,~,stat] = ttest2(corrvals1,corrvals2);
orig_tstat = stat.tstat;

alldata = cat(2,data1,data2);
origdesign = Make_designVect([size(data1,2) size(data2,2)]);


for i = 1:nperm
    design = origdesign(randperm(length(origdesign)));
    newcorrmat1 = corr(alldata(:,find(design == 1)),'Type',type);
    newcorrmat2 = corr(alldata(:,find(design == 2)),'Type',type);
    newcorrvals1 = newcorrmat1(find(belowDiag(ones(size(newcorrmat1)))));
    newcorrvals2 = newcorrmat2(find(belowDiag(ones(size(newcorrmat2)))));
    [~,~,~,stat] = ttest2(newcorrvals1,newcorrvals2);
    perm_tstat(i) = stat.tstat;
end

p = value_prctile(perm_tstat,orig_tstat);
