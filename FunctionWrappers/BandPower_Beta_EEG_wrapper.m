function [BetaPowerOut] = BandPower_Beta_EEG_handle(EEG)

disp(' ')
disp('Computing beta power...')
%EEG = pop_eegfiltnew(EEG, 7, 13, 826, 0, [], 0);
%[spectra,freqs] = spectopo(EEG.data,0,EEG.srate,'plot','off');
for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    BetaPowerOut(c) = trapz(EEG.spectra(c,intersect(find(EEG.freqs >= 13),find(EEG.freqs <= 30))));
end