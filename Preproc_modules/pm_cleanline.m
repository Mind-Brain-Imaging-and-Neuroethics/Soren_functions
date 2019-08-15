function preproc = pm_cleanline(cfg,data)

if isft(data)
   preproc = ft2eeglab(data); 
else
    preproc = data;
end

preproc = pop_cleanline(preproc, 'bandwidth',cfg.bandwidth,'chanlist',1:preproc.nbchan,'computepower',1,'linefreqs',cfg.linefreq ,'normSpectrum',1,'p',0.05,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
delete(gcf)

if isft(data)
   preproc = eeglab2fieldtrip(preproc,'preprocessing','none');
   preproc = restore_ft(data,preproc);
end