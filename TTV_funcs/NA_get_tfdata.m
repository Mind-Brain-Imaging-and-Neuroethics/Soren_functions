function [timefreq_data] = NA_get_tfdata(settings,data)

freqs = settings.tfparams.fbands;

timefreq_data = cell(1,length(freqs));
switch settings.tfparams.method
    case 'hilbert'
        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
        data = ft_resampledata(cfg,data);
        
        if isfield(settings.tfparams,'trials') && ~strcmpi(settings.tfparams.trials,'all')
            cfg = []; cfg.trials = settings.tfparams.trials;
            data = ft_preprocessing(cfg,data);
        end
        
        for q = 1:settings.nfreqs
            if ~isempty(freqs{q})
                cfg = [];
                if isnan(freqs{q}(1))
                    cfg.lpfilter = 'yes'; cfg.lpfreq = freqs{q}(2); cfg.lpinstabilityfix = 'split';
                else
                    cfg.bpfilter = 'yes'; cfg.bpfreq = freqs{q}; cfg.bpinstabilityfix = 'split';
                end
                timefreq_data{q} = ft_preprocessing(cfg,data);
            else
                timefreq_data{q} = data;
            end
            cfg = []; cfg.hilbert = 'complex';
            timefreq_data{q} = ft_preprocessing(cfg,timefreq_data{q});
        end
        
    case 'wavelet'
        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
        data = ft_resampledata(cfg,data);
        
        if isfield(settings.tfparams,'trials') && ~strcmpi(settings.tfparams.trials,'all')
            cfg = []; cfg.trials = settings.tfparams.trials;
            data = ft_preprocessing(cfg,data);
        end
        
        data_allrange = (settings.pseudo.prestim(1)-settings.srate/5):(settings.real.poststim(end));
        cfg = []; cfg.method = 'wavelet'; cfg.output = 'fourier'; cfg.foi = logspace(freqs{2}(1),freqs{end}(2),50);
        cfg.keeptrials = 'yes'; cfg.toi = data.time{1}(data_allrange); cfg.width = 3;
        freqdata = ft_freqanalysis(cfg,data);
        foi{i} = cfg.foi;
        
        timefreq_data{1} = data;
        for c = 1:length(timefreq_data{1}.trial)
            timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c}(:,data_allrange);
        end
        
        for c = 1:length(cfg.foi)
            for cc = 1:length(data.trial)
                timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,c,:));
            end
            timefreq_data{c+1}.time = freqdata.time;
            timefreq_data{c+1}.label = data.label;
            for cc = 1:length(freqs)
                if ~isempty(freqs{cc}) && cfg.foi(c) >= freqs{cc}(1) && cfg.foi(c) <= freqs{cc}(2)
                    timefreq_data{c+1}.parent = cc;
                end
            end
        end
        freqdata = [];
    case 'fft'
        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
        data = ft_resampledata(cfg,data);
        
        if isfield(settings.tfparams,'trials') && ~strcmpi(settings.tfparams.trials,'all')
            cfg = []; cfg.trials = settings.tfparams.trials;
            data = ft_preprocessing(cfg,data);
        end
        
        data_allrange = (settings.pseudo.prestim(1)-settings.srate/5):(settings.real.poststim(end));
        cfg = []; cfg.method = 'mtmconvol'; cfg.output = 'fourier'; cfg.foi = freqs{2}(1):2:freqs{end}(2);
        cfg.keeptrials = 'yes'; cfg.taper = 'hanning'; cfg.t_ftimwin = ones(length(cfg.foi))*1;
        cfg.toi = [data.time{1}(data_allrange)];
        freqdata = ft_freqanalysis(cfg,data);
        foi{i} = cfg.foi;
        
        timefreq_data{1} = data;
        for c = 1:length(timefreq_data{1}.trial)
            timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c}(:,data_allrange);
        end
        
        for c = 1:length(cfg.foi)
            for cc = 1:length(data.trial)
                timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,c,:));
            end
            timefreq_data{c+1}.time = freqdata.time;
            timefreq_data{c+1}.label = data.label;
            for cc = 1:length(freqs)
                if ~isempty(freqs{cc}) && cfg.foi(c) >= freqs{cc}(1) && cfg.foi(c) <= freqs{cc}(2)
                    timefreq_data{c+1}.parent = cc;
                end
            end
        end
        freqdata = [];
        
    case 'irasa'
        if isfield(settings.tfparams,'trials') && ~strcmpi(settings.tfparams.trials,'all')
            cfg = []; cfg.trials = settings.tfparams.trials;
            data = ft_preprocessing(cfg,data);
        end
        
        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
        newdata = ft_resampledata(cfg,data);
        
        data_allrange = (settings.pseudo.prestim(1)-settings.srate/5):(settings.real.poststim(end));
        
        cfg = []; cfg.oscifrac = settings.tfparams.oscifrac; cfg.winsize = 2;
        cfg.toi = newdata.time{1}(data_allrange); cfg.foi = freqs(2:end);
        freqdata = IRASA_tf(cfg,data);
        
        freqdata.fourierspctrm = freqdata.fourierspctrm + 0.00001j; % add a tiny imaginary component for compatibility
        
        timefreq_data{1} = data;
        for c = 1:length(timefreq_data{1}.trial)
            timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c}(:,data_allrange);
        end
        
        for c = 1:length(cfg.foi)
            for cc = 1:length(data.trial)
                timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,c,:));
            end
            timefreq_data{c+1}.time = freqdata.time;
            timefreq_data{c+1}.label = data.label;
        end
        freqdata = [];
        
end