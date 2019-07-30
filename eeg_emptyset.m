function EEG = eeg_emptyset

data = struct; 
data.trial{1} = zeros(1,1);
data.time{1} = 1;
data.fsample = 1;

EEG = ft2eeglab(data);