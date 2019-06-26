function [bp] = Bandpower_EEG_wrapper(EEG,frange,norm_bandpass)

if nargin < 3
   norm_bandpass = 'no'; 
end

for c = 1:EEG.nbchan
   bp(c) = bandpower(EEG.data(c,:),EEG.srate,frange); 
   if ~strcmpi(norm_bandpass,'no')
      bp(c) = bp(c)/bandpower(EEG.data(c,:),EEG.srate,norm_bandpass); 
   end
end