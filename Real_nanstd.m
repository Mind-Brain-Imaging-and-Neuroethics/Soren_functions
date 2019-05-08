function [outm] = Real_nanstd(inm,dimn)

if dimn == 2
    for c = 1:size(inm,1)
        outm(c) = nan_std(squeeze(inm(c,:)));
    end
elseif dimn == 3
    for c = 1:size(inm,1)
        for cc = 1:size(inm,2)
            outm(c,cc) = nan_std(squeeze(inm(c,cc,:)));
        end
    end
end
