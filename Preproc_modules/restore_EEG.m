function EEG = restore_EEG(origEEG,EEG)

fields = {'chanlocs','event','icaact','icawinv','icasphere','icaweights'};

for c = 1:length(fields)
   EEG.(fields{c}) = origEEG.(fields{c}); 
end

EEG = eeg_checkset(EEG);