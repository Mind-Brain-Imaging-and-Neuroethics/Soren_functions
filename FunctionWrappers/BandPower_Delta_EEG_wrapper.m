function [DeltaPowerOut] = BandPower_Delta_EEG_handle(EEG)


disp(' ')
disp('Computing delta power...')
%EEG = pop_eegfiltnew(EEG, 7, 13, 826, 0, [], 0);
%[spectra,freqs] = spectopo(EEG.data,0,EEG.srate,'plot','off');
for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    DeltaPowerOut(c) = trapz(EEG.spectra(c,intersect(find(EEG.freqs >= 1),find(EEG.freqs <= 4))));
end