function NA_analysis(settings)
%settings:
%  datasetname: string, name for the dataset
%  inputdir: string, directory where the relevant files are found - NO
%        slash at the end
%  outputdir: string, directory where calculated files are saved - NO slash
%        at the end
%  pseudo: struct with fields:
%        prestim: indices in the epoched data for prestim pseudotrial
%        poststim: indices in the epoched data for poststim pseudotrial
%  real: struct with fields:
%        prestim: indices in the epoched data for prestim real trial
%        poststim: indices in the epoched data for poststim real trial
%  aucindex: vector of data indices for computation of the TTV, ERSP, and
%        ITC indices
%  datasetinfo: a structure with the labels and channel/gradiometer
%        locations for the data set
%  datatype: 'EEG', 'MEG' or 'ECoG'
%  layout: .lay file or structure for plotting if MEG data, or EEGLAB
%        chanlocs structure for EEG data
%  steps: cell array containing which steps of the analysis to do
%        (default = {'all'})
%  srate: sample rate of the time-frequency data (default = 200)
%  units: 'prcchange', 'log', 'raw', or 'zscore' (default = 'prcchange')
%  files: the input to 'dir' to select what files to use (default = '*.mat'
%        for ECoG/MEG data, '*.set' for EEG data)
%  pool: the maximum parallel pool size (default = 24)
%  tfparams: a structure with fields relating to the frequency decomposition
%        method: 'hilbert', 'wavelet', or 'fft' - latter two not
%        implemented yet
%        fbands: cell array of frequency bands - use [] for broadband
%        fbandnames: names of the frequency bands of interest
%        condition: a number or numbers that indicate which condition to 
%        include. These numbers should be found in the data.trialinfo field
%        (only usable with fieldtrip data currently; default = 'all')
%        pf_adjust: adjust frequency bands for theta, alpha, and beta (if
%        input) using the subjects individual alpha frequency (default =
%        'yes')
%  rest: a structure with fields related to the resting state recordings
%        restdir: directory where resting state recordings are found
%        restfiles: the input to 'dir' to select what files to use (default
%        = '*.mat' for ECoG/MEG data, '*.set' for EEG data)
%        bandpass: a vector of the lowest and highest frequencies in the
%        data
%  behav: a cell array of structures containing the different behavioural
%        analyses to do. Each structure should have the fields:
%        file: a .mat file consisting of either a 1 by n subjects cell array
%        containing the behavioural variable of interest (implies a
%        within-subjects or mixed-effects design), or a n subjects by 1
%        vector. This should be either found in settings.outputdir, or the
%        whole path should be specified
%        indvar: either 'erspindex', 'ttvindex', 'rest_bp', or 'rest_relbp'
%        foi: can specify a frequency band name (default = 'all')


alltime_pre = cputime;
addpath('/group/northoff/share/fieldtrip-master')
ft_defaults
addpath('/group/northoff/share/eeglab14_1_2b')
addpath(genpath('/home/soren/Documents/MATLAB'))
rmpath(genpath('/home/soren/Documents/MATLAB/osl-core-master'))
rmpath(genpath('/home/soren/Documents/MATLAB/ImaGIN2'))

if strcmpi(settings.datatype,'EEG')
%    eeglab rebuild
%    addpath('/group/northoff/share/fieldtrip-master/external/eeglab')
end

settings = SetDefaults(settings);

cd(settings.inputdir)

if ~isempty(find(contains(settings.steps,'tf_filter')))
    disp('Performing time-frequency analysis...')
    tic;
    settings = NA_tf_func(settings);
    time_post = toc;
    disp(['Time-frequency analysis took ' num2str(time_post) ' seconds'])
end

cd(settings.outputdir)
% if ~isempty(find(contains(settings.steps,'calc')))
%     disp('Calculating TTV, ERSP, ITC...')
%     tic;
%     settings = TTV_ERSP_calc_func(settings);
%     time_post = toc;
%     disp(['Calculation of measures took ' num2str(time_post) ' seconds'])
% end

if ~isempty(find(contains(settings.steps,'results')))
    disp('Calculating indices and statistics...')
    tic;
    NA_results_func(settings)
    time_post = toc;
    disp(['Calculation took ' num2str(time_post) ' seconds'])
end

if ~isempty(find(contains(settings.steps,'restcorr')))
    disp('Calculating and correlating resting state parameters...')
    tic;
    NA_rest_func(settings)
    
    time_post = toc;
    disp(['Resting state analysis took ' num2str(time_post) ' seconds'])
end

if ~isempty(find(contains(settings.steps,'behav')))
    disp('Calculating behaviour relationships...')
    tic;
    NA_behaviour_func(settings)
    time_post = toc;
    disp(['Behaviour took ' num2str(time_post) ' seconds'])
end

if ~isempty(find(contains(settings.steps,'FDR')))
    disp('Correcting across analyses with FDR...')
    tic;
    NA_FDR_correct(settings)
    time_post = toc;
    disp(['FDR correction took ' num2str(time_post) ' seconds'])
end

if ~isempty(find(contains(settings.steps,'figures')))
    disp('Generating figures...')
    tic;
    system(['mkdir ' settings.outputdir '/' settings.datasetname '_figures'])
    NA_figures_func(settings)
    time_post = toc;
    disp(['Creating figures took ' num2str(time_post) ' seconds'])
end

alltime_post = cputime;

disp('Analysis finished')
disp(['All steps took ' num2str(alltime_post-alltime_pre) ' seconds'])
disp(['Results saved in ' settings.outputdir])

end


function settings = SetDefaults(settings)

%settings = settings;

settings = setdefault(settings,'comparefreqs','yes');

if ~isfield(settings,'files')
    if strmcpi(settings.datatype,'EEG')
        settings.files = '*.set';
    else
        settings.files = '*.mat';
    end
end

if isfield(settings,'rest')
    if strcmpi(settings.datatype,'EEG')
        settings.rest = setdefault(settings.rest,'restfiles','*.set');
    else
        settings.rest = setdefault(settings.rest,'restfiles','*.mat');
    end
end

settings.tfparams = setdefault(settings.tfparams,'pf_adjust','yes');

settings = setdefault(settings,'steps',{'all'});

if strcmpi(settings.steps{1},'all')
    settings.steps = {'tf_filter','results','restcorr','FDR','figures'};
end

settings = setdefault(settings,'srate',200);

settings = setdefault(settings,'units','prcchange');

settings.tfparams = setdefault(settings.tfparams,'trials','all');

settings.tfparams = setdefault(settings.tfparams,'method','hilbert');

settings.tfparams = setdefault(settings.tfparams,'continue','no');

settings = setdefault(settings,'fdr','yes');

settings = setdefault(settings,'load_allmeas','no');

if isfield(settings,'behav')
    for i = 1:length(settings.behav)
        settings.behav{i} = setdefault(settings.behav{i},'loadvar','no');
    end
end

settings.nfreqs = length(settings.tfparams.fbandnames);

settings = setdefault(settings,'pool',24);

if isfield(settings.datasetinfo,'label')
    settings.nbchan = length(settings.datasetinfo.label);
else
    settings.nbchan = length(settings.datasetinfo.atlas.tissuelabel);
end

end
