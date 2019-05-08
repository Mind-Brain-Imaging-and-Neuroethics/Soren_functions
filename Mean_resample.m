function [output] = Mean_resample(invector,sampfact)

for c = 1:length(invector)/sampfact
    output(c) = nan_mean(invector(((c-1)*sampfact+1):(c*sampfact)));
end