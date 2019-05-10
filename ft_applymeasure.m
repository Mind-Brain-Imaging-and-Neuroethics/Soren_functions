function [outputs] = ft_applymeasure(cfg)
% ft_applymeasure applies a certain measurement (ex. PLE, DFA, LZC) to a
% data set
%
% Input arguments:
%
% cfg: a configuration structure with the following fields and defaults
%
%      Basic options:
%
%      measure: a cell array containing the function handles to apply to
%         the data (required input). These function handles should take an
%         EEG structure, or at least a structure with fields 'data' (in the
%         format channels x samples) and 'srate' (the sampling rate).
%      format: 'eeglab' or 'fieldtrip' (default = 'eeglab');
%      dir: the base directory (default = get from user interface)
%      files: a string telling what files to select. Usually this
%         includes a wildcard character ('*'). For example, to select all
%         files in the folder ending in 'myfiles.mat', enter '*myfiles.mat'.
%         You can also set up nested directory structures this way. To
%         select all files named 'myfiles.mat' in all directories one level
%         down from the current one, enter */*myfiles.mat. (default =
%         '*.mat' for fieldtrip format, '*.set' for eeglab format)
%      outfile: the file into which to write the output. Can specify 'none'
%         to skip writing to disk (default = 'outputs.mat')
%
%      Data transformations:
%
%      filter: band-pass filters the data into the input range [low_cutoff,
%         high_cutoff] using a Butterworth filter (default = 'no')
%      hilbert: take the Hilbert envelope of the signal (required for DFA).
%         Uses the implementation in the NBT toolbox (default = 'no')
%      irasa: decompose the data using IRASA. By default uses a 10-second
%         window size, non-overlapping, with additional resampling factors
%         (as in Muthukumaraswamy and Liley, 2017). (default = 'no')
%
%      Other options:
%
%      continue: This option loads the output file and restarts the loop
%         from the last subject saved (default = 'no')
%      ftvar: for fieltrip format - the variable name in the .mat file
%         which corresponds to your data (default = the first structure
%         found in the data file)
%      concatenate: concatenate resting-states organized into trials
%         (default = 'yes')
%      parallel: a structure with the following options:
%         use_parallel: use parallel processing (default = 'no');
%         pool: size of parallel pool (default = system default pool size)
%      irasa_save_specs: for irasa - option to save spectra output from
%         IRASA to make future computations faster (default = 'no')
%
%      Surrogating and subsampling:
%
%      surrogate: a structure the following options
%         do_surr: do the surrogating or not (default = 'no');
%         nsurr: the number of surrogates to perform (default = 200);
%         method: surrogate method - choose either 'random', 'aaft', or
%         'iaaft' (default = 'aaft')
%      subsample: you can specify a sub-range of the data on which to apply
%         the measures. This is a structure with the following options
%
%         do_subsample: do the subsampling or not (default = 'no');
%         length: the length of the subsample (required input if you input
%            a subsample structure)
%         startpoint: a starting point from which to draw the
%            samples, or 'random' to randomly select one. (default =
%            'random')
%         bootstrap: sets the number of subsamples to take (default = 1);
%
%
% Outputs:
%
% outputs: a structure with the following fields:
%      data: the applied measures
%      dimord: a string telling you the order of the dimensions of 'data'.
%         For example, 'sub_chan_meas' means that the first dimension is
%         subjects, the second dimension is channels, and the third
%         dimension is measures.
%      sub: the file names that the function read
%      chan: the channel labels of the data files
%      meas: the measures you input in cfg.measure
%      elec or grad: the electrode or gradiometer position information -
%         used later in ft_measurestatistics



%% Set up defaults
if ~cfgcheck(cfg,'format')
    cfg.format = 'eeglab';
end

if ~cfgcheck(cfg,'files')
    if cfgcheck(cfg,'format','eeglab')
        cfg.files = '*.set';
    else
        cfg.files = '*.mat';
    end
end

if ~cfgcheck(cfg,'outfile')
    cfg.outfile = fullfile(pwd,'outputs.mat');
end

if ~cfgcheck(cfg,'concatenate')
    cfg.concatenate = 'yes';
end

if ~cfgcheck(cfg,'filter')
    cfg.filter = 'no';
end

if ~cfgcheck(cfg,'hilbert')
    cfg.hilbert = 'no';
end

if ~cfgcheck(cfg,'irasa')
    cfg.irasa = 'no';
end

if ~cfgcheck(cfg,'irasa_save_specs')
   cfg.irasa_save_specs = 'no'; 
end

if ~cfgcheck(cfg,'continue')
    cfg.continue = 'no';
end

if ~cfgcheck(cfg,'parallel')
    cfg.parallel.do_parallel = 'no';
end

if ~cfgcheck(cfg.parallel,'pool')
    cfg.parallel.pool = 'default';
end

if ~cfgcheck(cfg,'surrogate')
    cfg.surrogate.do_surr = 'no';
elseif cfgcheck(cfg,'surrogate','yes')
    if ~cfgcheck(cfg.surrogate,'nsurr')
        cfg.surrogate.nsurr = 200;
    end
    if ~cfgcheck(cfg.surrogate,'method')
        cfg.surrogate.method = 'aaft';
    end
end

if ~cfgcheck(cfg,'subsample')
    cfg.subsample.do_subsample = 'no';
elseif cfgcheck(cfg,'subsample','yes')
    if ~cfgcheck(cfg.subsample,'startpoint')
        cfg.subsample.startpoint = 'random';
    end
    if ~cfgcheck(cfg.subsample,'bootstrap')
        cfg.subsample.bootstrap = 1;
    end
end


%% Start up fieldtrip, eeglab if necessary
if cfgcheck(cfg,'format','eeglab')
    eeglab rebuild
end

ft_defaults

%% Find the base directory
dirname = cfgparse(cfg,'dir');
if isnan(dirname)
    dirname = uigetdir;
end


files = dir(fullfile(dirname,cfg.files));

%% Set up the outputs structure
if cfgcheck(cfg,'continue','yes')
    load(cfg.outfile)
else
    outputs = struct;
    if cfgcheck(cfg.subsample,'do_subsample','no') && cfgcheck(cfg.surrogate,'do_surr','no')
        outputs.dimord = 'sub_chan_meas';
    elseif cfgcheck(cfg.subsample,'do_subsample','no') && cfgcheck(cfg.surrogate,'do_surr','yes')
        outputs.dimord = 'sub_chan_meas_surr';
        outputs.surr = 1:cfg.surrogate.nsurr;
    elseif cfgcheck(cfg.subsample,'do_subsample','yes') && cfgcheck(cfg.surrogate,'do_surr','no')
        if cfg.subsample.bootstrap > 1
            outputs.dimord = 'sub_chan_meas_sample';
            outputs.sample = 1:cfg.subsample.bootstrap;
        else
            outputs.dimord = 'sub_chan_meas';
        end
        %     elseif cfgcheck(cfg.subsample,'do_subsample','yes') && cfgcheck(cfg.surrogate,'do_surr','yes')
        %         outputs.dimord = 'sub_chan_meas_sample_surr';
    end
    outputs.meas = cfg.measure;
    outputs.startsub = 1;
    
end

if cfgcheck(cfg.parallel,'use_parallel','no')
    
    for i = outputs.startsub:length(files)
        
        %% Load the data
        filename = files(i).name;
        outputs.sub{i} = files(i).name;
        outputs.startsub = i+1;
        
        disp(' ')
        disp(['Now processing subject ' num2str(i)])
        
        
        if cfgcheck(cfg,'format','fieldtrip')
            EEG = struct;
            allvars = load(fullfile(files(i).folder,filename));
            if ~cfgcheck(cfg,'ftvar')
                names = fieldnames(allvars);
                for c = 1:length(names)
                    if isstruct(allvars.(names{c}))
                        data = allvars.(names{c});
                        clear allvars
                        break
                    end
                end
                clear names
            else
                data = allvars.(ftvar);
            end
            
            if cfgcheck(cfg,'concatenate','yes')
                data = ft_concat(data);
            end
            
            if i == outputs.startsub-1
                if isfield(data,'label')
                    outputs.chan = data.label;
                end
                if isfield(data,'elec')
                    outputs.elec = data.elec;
                elseif isfield(data,'grad')
                    outputs.grad = data.grad;
                end
            end
            
            EEG = ft2eeglab(data);
        else
            EEG = pop_loadset( 'filename', filename, 'filepath', files(i).folder);
            outputs.chanlocs = EEG.chanlocs;
            if i == outputs.startsub-1
                data = eeglab2fieldtrip(EEG,'preprocessing','none');
                if isfield(data,'label')
                    outputs.chan = data.label;
                end
                if isfield(data,'elec')
                    outputs.elec = data.elec;
                elseif isfield(data,'grad')
                    outputs.grad = data.grad;
                end
            end
        end
        
        EEG.filename = fullfile(files(i).folder,files(i).name);

        
        %% Subsample, surrogate, and apply the measures
        if cfgcheck(cfg.surrogate,'do_surr','no') && cfgcheck(cfg.subsample,'do_subsample','no')
            EEG = transform_data(cfg,EEG);
            for c = 1:length(cfg.measure)
                outputs.data(i,:,c) = cfg.measure{c}(EEG);
            end
            
        elseif cfgcheck(cfg.subsample,'do_subsample','yes') && cfgcheck(cfg.surrogate,'do_surr','no')
            disp(['Taking ' num2str(cfg.subsample.bootstrap) ' subsamples...'])
            for q = 1:cfg.subsample.bootstrap
                                tmpcfg = cfg;
                tmpcfg.irasa = 'no';
                EEG = transform_data(cfg,EEG);
                
                newEEG = SubSample(cfg,EEG);
                
                if cfgcheck(cfg,'irasa','yes')
                    tmpcfg = cfg;
                    tmpcfg.hilbert = 'no';
                    tmpcfg.filter = 'no';
                    EEG = transform_data(tmpcfg,EEG);
                end
                
                fprintf([num2str(q) ' '])
                for c = 1:length(cfg.measure)
                    outputs.data(i,:,q,c) = cfg.measure{c}(newEEG);
                end
            end
        elseif cfgcheck(cfg.surrogate,'do_surr','yes') && cfgcheck(cfg.subsample,'do_subsample','no')
            for c = 1:EEG.nbchan
                disp(' ')
                disp(['Creating surrogates for channel ' num2str(c)])
                if ~cfgcheck(cfg,'filter','no')
                    tmpcfg = cfg; 
                    tmpcfg.irasa = 'no';
                    tmpcfg.hilbert = 'no';
                    EEG = transform_data(tmpcfg,EEG);
                end
                
                switch cfg.surrogate.method
                    case 'aaft'
                        tmp = AAFT(EEG.data(c,:),cfg.surrogate.nsurr);
                    case 'iaaft'
                        tmp = IAAFT(EEG.data(c,:),cfg.surrogate.nsurr);
                    case 'random'
                        tmp = zeros(EasyParse(varargin,'Surrogate'),length(EEG.data));
                        for q = 1:cfg.surrogate.nsurr
                            tmp(q,:) = EEG.data(c,randperm(length(EEG.data)));
                        end
                end
                newEEG = EEG;
                newEEG.data = tmp';
                newEEG.nbchan = cfg.surrogate.nsurr;
                
                if cfgcheck(cfg,'irasa','yes') || cfgcheck(cfg,'hilbert','yes')
                    tmpcfg = cfg; 
                    tmpcfg.filter = 'no';
                    EEG = transform_data(tmpcfg,EEG);
                end
                
                for cc = 1:length(cfg.measure)
                    outputs.data(i,c,cc,:) = cfg.measure{cc}(newEEG);
                end
            end
        elseif cfgcheck(cfg.surrogate,'do_surr','yes') && cfgcheck(cfg.subsample,'do_subsample','yes')
            % finish later
        end
        
        %% Save after every subject so you can continue later
        if ~cfgcheck(cfg,'outfile','none')
            try
                save(cfg.outfile,'outputs');
            catch
                warning('Saving failed')
            end
        end
        
    end
    
else
    %% Parallel version
    
    if cfgcheck(cfg,'format','fieldtrip')
        allvars = parload(fullfile(files(1).folder,filename));
        if ~cfgcheck(cfg,'ftvar')
            names = fieldnames(allvars);
            for c = 1:length(names)
                if isstruct(allvars.(names{c}))
                    tmpdata = allvars.(names{c});
                    allvars = [];
                    break
                end
            end
            names = [];
        else
            tmpdata = allvars.(ftvar);
        end
    else
        EEG = pop_loadset('filename',files(1).name,'filepath',files(1).folder);
        outputs.chanlocs = EEG.chanlocs;
        tmpdata = eeglab2fieldtrip(EEG,'preprocessing','none');
    end
    
    if isfield(tmpdata,'label')
        outputs.chan = tmpdata.label;
    end
    if isfield(tmpdata,'elec')
        outputs.elec = tmpdata.elec;
    elseif isfield(tmpdata,'grad')
        outputs.grad = tmpdata.grad;
    end
    
    if ~cfgcheck(cfg.parallel,'pool','default')
        delete(gcp('nocreate'))
        parpool(cfg.parallel.pool)
    end
    
    clear data
    
    sub = cell(1,length(files));
    outdata = cell(1,length(files));
    
    
    parfor i = outputs.startsub:length(files)
        
        %% Load the data
        filename = files(i).name;
        sub{i} = files(i).name;
        %outputs.startsub = i+1;
        
        disp(' ')
        disp(['Now processing subject ' num2str(i)])
        
        
        if cfgcheck(cfg,'format','fieldtrip')
            EEG = struct;
            allvars = parload(fullfile(files(i).folder,filename));
            if ~cfgcheck(cfg,'ftvar')
                names = fieldnames(allvars);
                for c = 1:length(names)
                    if isstruct(allvars.(names{c}))
                        data = allvars.(names{c});
                        allvars = [];
                        break
                    end
                end
                names = [];
            else
                data = allvars.(ftvar);
            end
            
            if cfgcheck(cfg,'concatenate','yes')
                data = ft_concat(data);
            end
            
            EEG = ft2eeglab(data);
        else
            EEG = pop_loadset( 'filename', filename, 'filepath', files(i).folder);
        end
        
        EEG.filename = fullfile(files(i).folder,files(i).name);
        
        %% Subsample, surrogate, and apply the measures
        if cfgcheck(cfg.surrogate,'do_surr','no') && cfgcheck(cfg.subsample,'do_subsample','no')
            EEG = transform_data(cfg,EEG);
            for c = 1:length(cfg.measure)
                outdata{i}(1,:,c) = cfg.measure{c}(EEG);
            end
            
        elseif cfgcheck(cfg.subsample,'do_subsample','yes') && cfgcheck(cfg.surrogate,'do_surr','no')
            disp(['Taking ' num2str(cfg.subsample.bootstrap) ' subsamples...'])
            for q = 1:cfg.subsample.bootstrap
                tmpcfg = cfg;
                tmpcfg.irasa = 'no';
                EEG = transform_data(cfg,EEG);
                
                newEEG = SubSample(cfg,EEG);
                
                if cfgcheck(cfg,'irasa','yes')
                    tmpcfg = cfg;
                    tmpcfg.hilbert = 'no';
                    tmpcfg.filter = 'no';
                    EEG = transform_data(tmpcfg,EEG);
                end
                
                fprintf([num2str(q) ' '])
                EEG = transform_data(cfg,EEG);
                for c = 1:length(cfg.measure)
                    outdata{i}(1,:,q,c) = cfg.measure{c}(newEEG);
                end
            end
        elseif cfgcheck(cfg.surrogate,'do_surr','yes') && cfgcheck(cfg.subsample,'do_subsample','no')
            for c = 1:EEG.nbchan
                disp(' ')
                disp(['Creating surrogates for channel ' num2str(c)])
                if ~cfgcheck(cfg,'filter','no')
                    tmpcfg = cfg; 
                    tmpcfg.irasa = 'no';
                    tmpcfg.hilbert = 'no';
                    EEG = transform_data(tmpcfg,EEG);
                end
                
                switch cfg.surrogate.method
                    case 'aaft'
                        tmp = AAFT(EEG.data(c,:),cfg.surrogate.nsurr);
                    case 'iaaft'
                        tmp = IAAFT(EEG.data(c,:),cfg.surrogate.nsurr);
                    case 'random'
                        tmp = zeros(EasyParse(varargin,'Surrogate'),length(EEG.data));
                        for q = 1:cfg.surrogate.nsurr
                            tmp(q,:) = EEG.data(c,randperm(length(EEG.data)));
                        end
                end
                newEEG = EEG;
                newEEG.data = tmp';
                newEEG.nbchan = cfg.surrogate.nsurr;
                
                if cfgcheck(cfg,'irasa','yes') || cfgcheck(cfg,'hilbert','yes')
                    tmpcfg = cfg; 
                    tmpcfg.filter = 'no';
                    EEG = transform_data(tmpcfg,EEG);
                end
                
                for cc = 1:length(cfg.measure)
                    outdata{i}(1,c,cc,:) = cfg.measure{cc}(newEEG);
                end
            end
        elseif cfgcheck(cfg.surrogate,'do_surr','yes') && cfgcheck(cfg.subsample,'do_subsample','yes')
            % finish later
        end
        
        
    end
    
    outputs.data = cat(1,outdata{:});
    outputs.sub = sub;
    
    if ~cfgcheck(cfg,'outfile','none')
        try
            parsave(cfg.outfile,'outputs',outputs);
        catch
            warning('Saving failed')
        end
    end
    
end


end

%%
function [EEG] = SubSample(cfg,EEG)
sampleSize = cfg.subsample.length;

if cfgcheck(cfg.subsample,'startpoint','random') %if no startpoint specified, pick a random one
    startPoint = randi(length(EEG.data)-sampleSize);
else %specify a specific latency to start at
    startPoint = cfg.subsample.startpoint;
end

disp(' ')
disp(['Processing data from data point ' num2str(startPoint) ' to data point ' num2str(startPoint+sampleSize) '...'])

if startPoint+sampleSize < length(EEG.data)
    EEG.data = EEG.data(:,startPoint:(startPoint+sampleSize));
else
    warning('Not enough samples to meet sample size requirement')
    disp(' ')
    disp(['Continuing with ' num2str(length(EEG.data)-startPoint) ' samples...'])
    disp([num2str(startPoint+sampleSize-length(EEG.data)) ' samples missing...'])
    EEG.data = EEG.data(:,startPoint:end);
end
end

function [EEG] = transform_data(cfg,EEG)

%% Apply transformation options

if ~cfgcheck(cfg,'filter','no')
    data = eeglab2fieldtrip(EEG,'preprocessing','none');
    
    tmpcfg = []; tmpcfg.bpfilter = 'yes'; tmpcfg.bpfreq = cfg.filter;
    tmpcfg.bpinstabilityfix = 'split';
    
    data = ft_preprocessing(tmpcfg,data);
    
    EEG = ft2eeglab(data);
end

if cfgcheck(cfg,'amplitude','yes')
    disp(' ')
    disp('Getting amplitude envelope...')
    Signal = EEG.data;
    
    Signal = transpose(Signal); %now channels are columns, time is rows. Needed for the DFA function
    
    SignalInfo = nbt_Info; %this initializes an Info Object
    SignalInfo.converted_sample_frequency = EEG.srate;
    
    bpfreqs = EasyParse(varargin,'AmpEn');
    AmplitudeEnvelope = nbt_GetAmplitudeEnvelope(Signal, SignalInfo, bpfreqs(1), bpfreqs(2), 2*(1/bpfreqs(1)));
    
    EEG.data = AmplitudeEnvelope';
end

if cfgcheck(cfg,'irasa','yes')
    disp(' ')
    disp('Performing IRASA...')
    if exist([EEG.filename '_IRASA_specs.mat'],'file')
        EEG = parload(EEG.filename,'EEG');
    else
        EEG = IRASA_window(EEG.data,EEG.srate);
    end
    
    if cfgcheck(cfg,'irasa_save_specs','yes')
       save([EEG.filename '_IRASA_specs.mat'],'EEG');
    end
end
end