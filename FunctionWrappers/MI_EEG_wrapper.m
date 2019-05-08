function [MI] = MI_EEG_wrapper(EEG)

for c = 1:size(EEG.data,1)
    tmp = mi(double(EEG.data(c,1:100)),double(EEG.data(c,(end-99):end)),'silent');
    MI(c) = tmp(2,1);
end
