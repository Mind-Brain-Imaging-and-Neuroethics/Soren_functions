function [bp] = Rel_Bandpower_EEG_wrapper(EEG,frange)

fields = fieldnames(EEG);
if contains(fields,'nbchan')
    for c = 1:EEG.nbchan
        bp(c) = bandpower(EEG.data(c,:),EEG.srate,frange)/bandpower(EEG.data(c,:),EEG.srate,[0.5 50]);
    end
else
    for c = 1:size(EEG.osci,2)
       bp(c) = bandpower(EEG.osci(:,c),EEG.freq,frange,'psd')/bandpower(EEG.osci(c,:),EEG.freq,[0.5 50],'psd');
    end
end
    