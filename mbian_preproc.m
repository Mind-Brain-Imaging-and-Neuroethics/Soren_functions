function [data] = mbian_preproc(cfg,data)
% mbian_preproc implements a standard preprocessing pipeline for the Mind,
% Brain Imaging and Neuroethics group at the Royal
%
% Input arguments:
%      cfg: a config structure with the following fields
%           datatype: 'eeg' or 'meg' (default = 'eeg')
%           format: 'eeglab' or 'fieldtrip' (default = 'eeglab' for eeg,
%              'fieldtrip for meg)
%           bandpass: [low_freq high_freq] (default = [0.5 50] for eeg,
%              [0.05 150] for meg
%           resample: number specifying the sample rate for the
%              preprocessed recording (default = 500 Hz or the original
%              sampling rate, whichever is lower)
%
%           Datatype-specific options:
%
%           EEG:
%
%           chanlocs: either 'lookup' to look up
%              channel locations in standard file, a string specifying the
%              file to look up locations in, or an EEGLAB chanlocs
%              structure (default = 'lookup')
%           line: structure with the following fields:
%              method: Can be 'cleanline', 'none', or 'notch' (default =
%                 'cleanline')
%              freq: an n x m matrix of line noise frequencies. For method
%                 'cleanline', input only the actual line frequencies (m =
%                 1). For method 'notch', input bands you want to notch
%                 filter (m = 2).
%           asr: clean data using artifact subspace reconstruction - 'yes'
%              or 'no' (default = 'yes')
%           reference: referencing scheme - 'default','average', or 'REST'
%              (default = 'average')
%           ica: decompose data using ICA - 'yes' or 'no' (default = 'yes')
%           mara: automatic ICA component rejection using MARA - 'yes' or
%              'no' (default = 'yes')
%
%           MEG:
%
%           not implemented for now
%
%      data: an EEGLAB or Fieldtrip data structure
%
% Outputs:
%      data: the preprocessed data structure
%
%

%% Set up defaults
cfg = setdefault(cfg,'datatype','eeg');
if ~cfgcheck(cfg,'format')
    if cfgcheck(cfg,'datatype','eeg')
        cfg.format = 'eeglab';
    elseif cfgcheck(cfg,'datatype','meg')
        cfg.format = 'fieldtrip';
    end
end

if cfgcheck(cfg,'datatype','eeg')
    cfg = setdefault(cfg,'bandpass',[0.5 50]);
elseif cfgcheck(cfg,'datatype','meg')
    cfg = setdefault(cfg,'bandpass',[0.05 150]);
end

if cfgcheck(cfg,'datatype','eeg') && data.srate > 500
    cfg = setdefault(cfg,'resample',500);
elseif cfgcheck(cfg,'datatype','meg') && data.fsample > 500
    cfg = setdefault(cfg,'resample',500);
end

if cfgcheck(cfg,'datatype','eeg')
    cfg = setdefault(cfg,'chanlocs','lookup');
    cfg = setdefault(cfg,'line','cleanline');
    cfg = setdefault(cfg,'asr','yes');
    cfg = setdefault(cfg,'reference','average');
    cfg = setdefault(cfg,'ica','yes');
    cfg = setdefault(cfg,'mara','yes');
elseif cfgcheck(cfg,'datatype','meg')
    % add later
end

ft_defaults

%% EEG pipeline
if cfgcheck(cfg,'datatype','eeg')
    eeglab rebuild
    
    % Read in channel locations
    if isstr(cfg.chanlocs)
        if cfgcheck(cfg,'chanlocs','lookup')
            EEG = pop_chanedit(EEG,'lookup',fullfile('plugins','dipfit2.3','standard_BESA','standard-10-5-cap385.elp'),'eval','chans = pop_chancenter( chans, [],[]);');
        else
            EEG = pop_chanedit(EEG,'lookup',cfg.chanlocs,'eval','chans = pop_chancenter( chans, [],[]);');
        end
    else
        EEG.chanlocs = cfg.chanlocs;
    end
    
    EEG  = eeg_checkset(EEG);
    
    
    % Filter the data
    % filter in fieldtrip so you don't have to manually set the filter order
    chanlocs = EEG.chanlocs; %ft2eeglab doesn't handle chanlocs well, so save these
    
    data = eeglab2fieldtrip(EEG,'preprocessing','none');
    tmpcfg = [] ; tmpcfg.bpfilter = 'yes'; tmpcfg.bpfreq = cfg.bandpass; tmpcfg.bpfilttype = 'fir';
    data = ft_preprocessing(tmpcfg,data);
    EEG = ft2eeglab(data);
    EEG.chanlocs = chanlocs;
    
    EEG  = eeg_checkset(EEG);
    
    
    % Clean line noise
    switch cfg.line.method
        case 'cleanline'
            EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',1:EEG.nbchan,'computepower',1,'linefreqs',cfg.line.freq ,'normSpectrum',1,'p',0.05,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
            delete(gcf)
        case 'notch'
            chanlocs = EEG.chanlocs; %ft2eeglab doesn't handle chanlocs well, so save these
            
            data = eeglab2fieldtrip(EEG,'preprocessing','none');
            tmpcfg = [] ; tmpcfg.bsfilter = 'yes'; tmpcfg.bsfreq = cfg.line.freq; tmpcfg.bpfilttype = 'fir';
            data = ft_preprocessing(tmpcfg,data);
            EEG = ft2eeglab(data);
            EEG.chanlocs = chanlocs;
    end
    EEG  = eeg_checkset(EEG);
    
    
    % Apply ASR
    if cfgcheck(cfg,'asr','yes')
        orig_chanlocs = EEG.chanlocs;
        
        EEG = clean_rawdata(EEG, 5, [-1], 0.8, -1, 5, 0.5);
        
        EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
        
        EEG  = eeg_checkset(EEG);
        
    end
    
    % Rereference
    switch cfg.reference
        case 'average'
            EEG = pop_reref( EEG, []);
        case 'REST'
            % add REST reference
    end
    EEG  = eeg_checkset(EEG);
    
    
    % Run ICA
    if cfgcheck(cfg,'ica','yes')
        EEG = pop_runica(EEG, 'interupt','off');
        EEG  = eeg_checkset(EEG);
        
    end
    
    % Reject components with MARA
    if cfgcheck(cfg,'mara','yes')
        artcomps = MARA(EEG);
        EEG  = pop_subcomp(EEG,artcomps,0);
        EEG  = eeg_checkset(EEG);
    end
    
elseif cfgcheck(cfg,'datatype','meg')
    %% MEG pipeline
end


