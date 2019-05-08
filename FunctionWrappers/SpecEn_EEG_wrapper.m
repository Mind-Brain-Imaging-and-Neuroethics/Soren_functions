function [SpecEnOut] = SpecEn_EEG_handle(EEG)

SpecEnOut = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing spectral entropy...')

for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    SpecEnOut(c) = pentropy(EEG.data(c,:),EEG.srate,'Instantaneous',false);
end