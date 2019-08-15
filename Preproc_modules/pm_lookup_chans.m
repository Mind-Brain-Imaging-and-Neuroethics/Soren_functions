function preproc = pm_lookup_chans(cfg,data)

if isft(data)
    preproc = ft2eeglab(data);
    preproc = ft_eeg_setlabels(preproc,EEG);
else
    preproc = data;
end

eegdir = extractBefore(which('eeglab'),'eeglab.m');

preproc = pop_chanedit(preproc,'lookup',fullfile(char(eegdir),'plugins','dipfit2.3','standard_BESA','standard-10-5-cap385.elp'),'eval','chans = pop_chancenter( chans, [],[]);');

if isft(data)
   preproc = eeglab2fieldtrip(preproc,'preprocessing','none'); 
   preproc = restore_ft(data,preproc);
end