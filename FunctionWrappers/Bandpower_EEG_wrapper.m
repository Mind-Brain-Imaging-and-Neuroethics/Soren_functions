function [bp] = Bandpower_EEG_wrapper(EEG,frange)

for c = 1:EEG.nbchan
   bp(c) = bandpower(EEG.data(c,:),EEG.srate,frange); 
end