function preproc = pm_averef(cfg,data)

if isEEG(data)
   preproc = pop_reref(data,[]); 
else
   cfg.reref = 'yes'; cfg.refchannel = 'all'; cfg.refmethod = 'avg';
   preproc = ft_preprocessing(cfg,data);
end