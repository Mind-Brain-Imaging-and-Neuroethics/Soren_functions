function [data,cont_data] = camcan_preproc(subid,filename,cont_data)

% Assuming ComputeCanada
basedir = extractBefore(filename,['sub-' subid]);
%basedir = '/scratch/sorenwt/camcan/cc700/mri/pipeline/release004/BIDSsep/megmax_task';

% Option to input the continuous data in case something goes wrong and you
% only want to do the epoched preprocessing
if ~exist('cont_data','var')
    
    %% First do movement compensation in MNE somehow
    
    %% Load in the file
    rawfile = filename;
    hdr = ft_read_header(rawfile);
    
    
    cfg = []; cfg.dataset = rawfile;
    cfg.hdr = rawfile;
    cont_data = ft_preprocessing(cfg);
    
    %% Filter
    
    cfg = []; cfg.lpfilter = 'yes'; cfg.lpfreq = 200;
    cfg.bsfilter = 'yes'; cfg.bsfreq = [49 51; 99 101; 149 151; 199 201];
    cont_data = ft_preprocessing(cfg,cont_data);
    
    %% Use autoreject to find threshold to exclude from ICA
    cfg = []; cfg.event = 1:2*cont_data.fsample:floor(length(cont_data.time{1}));
    cfg.event(end) = []; cfg.epoch = [0 (2*cont_data.fsample)-1];
    data = ft_epoch(cfg,cont_data);
    data.trialinfo = ones(length(data.sampleinfo),1);
    
    save(fullfile(basedir,['sub-' subid],[subid '_cont_epochs.mat']),'data')
    
    pyscript = fopen([subid '_pyscript.py'],'w');
    fprintf(pyscript,'import sys \n')
    fprintf(pyscript,"sys.path.insert(0, '/home/sorenwt/projects/def-gnorthof/sorenwt/MATLAB/Functions') \n")
    fprintf(pyscript,'from mne_preproc import autoreject_log \n')
    fprintf(pyscript,['autoreject_log(''' rawfile, ''',''' fullfile(basedir,['sub-' subid],[subid '_cont_epochs.mat'])...
        ''',''' fullfile(basedir,['sub-' subid],[subid '_badsegs.json']) ''')'])
    system(['python -c ''' subid '_pyscript.py'''])
    system(['rm ' subid '_pyscript.py'])
    
    bads = jsonread(fullfile(basedir,['sub-' subid],[subid '_badsegs.json']));

    %% Remove bad segments and re-concatenate, filter at 1 Hz for ICA
    
    cfg = []; cfg.trials = ~bads;
    data = ft_selectdata(cfg,data);
    
    cont_data_clean = data;
    cont_data_clean.trial = []; 
    cont_data_clean.trial{1} = cat(2,data.trial{:});
    cont_data_clean.time = []; 
    cont_data_clean.time{1} = linspace(0,length(cont_data_clean.trial{1})/data.fsample,length(cont_data_clean.trial{1}));
    data = [];
    
    cfg = []; cfg.hpfilter = 'yes'; cfg.hpfreq = 1; 
    cont_data_clean = ft_preprocessing(cfg,cont_data_clean);
    
    %% ICA using OSL-AFRICA
    
    D = spm_eeg_ft2spm(cont_data_clean,fullfile(basedir,['sub-' subid],'tmp'));
    D = osl_africa(D,'used_maxfilter',true); 
    system(['rm ' fullfile(basedir,['sub-' subid],'tmp*')])
%     cont_data_clean = spm2fieldtrip(D);
%     cont_data_clean.hdr = hdr; cont_data_clean.fsample = hdr.Fs;
    cont_data_clean = [];

    %% Apply ICA matrix back to original data
    
    cfg = []; cfg.unmixing = D.ica.sm'; cfg.topolabel = data.label;
    comp = ft_componentanalysis(cfg,cont_data);
    
    cfg = []; cfg.component = D.ica.bad_components; 
    cont_data = ft_rejectcomponent(cfg,cont_data);
    
end

%% Epoch data around stimuli

event = ft_read_event(rawfile);
types = extractfield(event,'type');
tpoints = extractfield(event,'sample');
latencies = tpoints(find(strcmpi(types,'Trigger')));
allvalues = extractfield(event,'value');
values = allvalues(find(strcmpi(types,'Trigger')));
latencies(find(values < 0)) = []; 

cfg = []; cfg.event = latencies; cfg.epoch = [-1.5*cont_data.fsample 1*cont_data.fsample]; %Epochs large in order to have trial padding
data = ft_epoch(cfg,cont_data);
data.trialinfo = values(find(values > 0));

%% Use autoreject to remove or interpolate bad epochs
save(fullfile(basedir,['sub-' subid],[subid '_ftdata.mat']),'data');
data = [];

pyscript = fopen([subid '_pyscript.py'],'w');
fprintf(pyscript,'import sys \n')
fprintf(pyscript,"sys.path.insert(0, '/home/sorenwt/projects/def-gnorthof/sorenwt/MATLAB/Functions') \n")
fprintf(pyscript,'from mne_preproc import autoreject_epochs \n')
fprintf(pyscript,['autoreject_epochs(''' rawfile, ''',''' fullfile(basedir,['sub-' subid],[subid 'cont_epochs_.mat'])...
    ''',''' fullfile(basedir,['sub-' subid],[subid '_mne.fif']) ''')'])
system(['python -c ' subid '_pyscript.py'])
system(['rm ' subid '_pyscript.py'])
    

cfg = []; cfg.dataset = fullfile(basedir,['sub-' subid],[subid '_mne.fif']);
cfg.hdr = cfg.dataset;
data = ft_preprocessing(cfg);
