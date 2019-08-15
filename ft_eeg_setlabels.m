function EEG = ft_eeg_setlabels(data,EEG)

for c = 1:length(data.label)
   EEG.chanlocs(c).labels = data.label{c}; 
end