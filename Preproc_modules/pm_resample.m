function preproc = pm_resample(cfg,data)

if isEEG(data)
   preproc = pop_resample(EEG,cfg.resamplefs);
else
    preproc = ft_resampledata(cfg,data);
end