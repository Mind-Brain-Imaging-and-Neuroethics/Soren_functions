function [pvals] = mcompare_bfholm(datain,statfun)

for c = 1:size(datain,2)
   for cc = 1:size(datain,2)
      pvals(c,cc) = statfun(datain(:,c),datain(:,cc)); 
   end
end


pvals = bonf_holm(pvals);

pvals(find(pvals > 1)) = 1;