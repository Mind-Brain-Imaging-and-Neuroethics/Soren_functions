function [pvalue] = bootdiff(funcn,samp1,samp2,nboot)

newfunc = @(x1,x2) funcn(x1)/funcn(x2); % difference is significant if the ratio is significantly different from 1

bootstat = zeros(1,nboot);

sampstat = newfunc(samp1,samp2);

for c = 1:nboot
   resamp1 = samp1(ceil(rand(size(samp1))*length(samp1)));
   resamp2 = samp2(ceil(rand(size(samp2))*length(samp2)));
   bootstat(c) = newfunc(resamp1,resamp2);
end

if sampstat > 1
   pvalue = length(find(bootstat <= 1))/length(bootstat)*2; 
else
    pvalue = length(find(bootstat >= 1))/length(bootstat)*2;
end

if pvalue > 1
    pvalue = 1;
end