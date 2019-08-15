function preproc = pm_asr(cfg,data)

if isft(data)
    preproc = ft2eeglab(data);
else
    preproc = data;
end

orig_chanlocs = preproc.chanlocs;

preproc = clean_rawdata(preproc, 5, [-1], 0.85, 4, 20,-1,'availableRAM_GB',4);

if isEEG(data)
    preproc = pop_interp(preproc, orig_chanlocs, 'spherical');
else
    preproc = eeglab2fieldtrip(EEG,'preprocessing','none');
    
    labels = preproc.label;
    
    preproc = restore_ft(data,preproc);
    
    if length(preproc.label) < length(data.label)
        tmpcfg = []; tmpcfg.method = 'spline'; tmpcfg.missingchannel = setdiff(preproc.label,data.label);
        tmpcfg2.method = 'triangulation'; tmpcfg.neighbours = ft_prepare_neighbours(tmpcfg2,preproc);
        preproc = ft_channelrepair(tmpcfg,preproc);
    end
end
