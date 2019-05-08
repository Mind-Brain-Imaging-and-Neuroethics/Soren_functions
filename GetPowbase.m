function [mean_powbase] = GetPowbase(EEG,windowLen,electrode)
%Inputs: EEG - the whole struct. 
%Windowlen - the length of the epoch you want a base for, in samples
%
%Outputs: mean_powbase

datalen = length(EEG.data);

mean_powbase = [];
for c = 1:floor(EEG.xmax/(windowLen/EEG.srate))
    [~,~,powbaserest{c}] = newtimef(EEG.data(electrode,1:datalen-mod(datalen,windowLen)),windowLen,[(c-1)*(1000*windowLen/EEG.srate)+1 c*(1000*windowLen/EEG.srate)],EEG.srate,0,'plotitc','off','plotersp','off');
    mean_powbase = [mean_powbase;powbaserest{c}];
end
mean_powbase = mean(mean_powbase,1);