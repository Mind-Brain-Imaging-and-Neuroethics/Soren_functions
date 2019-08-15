function preproc = pm_filter(cfg,data)

if isEEG(data)
   preproc = eeglab2fieldtrip(data,'preprocessing','none');
else
    preproc = data;
end

preproc = ft_preprocessing(cfg,preproc);

if isEEG(data)
   preproc = ft2eeglab(preproc);
   preproc = restore_EEG(data,preproc);
end

