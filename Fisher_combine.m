function [ptot] = Fisher_combine(pvals,varargin)

if CheckInput(varargin,'MCorrect')
    if EasyParse(varargin,'MCorrect','fdr')
        pvals = mafdr(pvals,'BHFDR',true);
    elseif EasyParse(varargin,'MCorrect','bonferroni')
        pvals = pvals*length(pvals);
    end
end

chi2 = -2*sum(log(pvals));

ptot = 1 - chi2cdf(chi2,length(pvals)*2);