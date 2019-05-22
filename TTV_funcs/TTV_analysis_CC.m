function TTV_analysis_CC(settings)
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
%        locations for the data set (default = extracted from preprocessing)
%  datatype: 'EEG', 'MEG' or 'ECoG'
%  layout: .lay file or structure for plotting if MEG data, or EEGLAB
%        chanlocs structure for EEG data
%  steps: cell array containing which steps of the analysis to do 
%        (default = {'all'})
%  srate: sample rate of the time-frequency data (default = 200)
%  units: 'prcchange', 'log', 'raw', or 'zscore' (default = 'prcchange')
%  tfparams: a structure with fields relating to the frequency decomposition
%        method: 'hilbert', 'wavelet', or 'fft' - latter two not
%        implemented yet
%        fbands: cell array of frequency bands - use [] for broadband
%        fbandnames: names of the frequency bands of interest
%        trials: indices of the trials of interest (default = 'all');
%  rest: a structure with fields related to the resting state recordings
%        restdir: directory where resting state recordings are found
%        preprocparams: a structure with fields related to the
%        preprocessing of the resting states. If this field is not
%        specified, it is assumed the files are already preprocessed in
%        eeglab format
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
addpath('/project/def-gnorthof/sorenwt/fieldtrip-master')
ft_defaults
addpath('/project/def-gnorthof/sorenwt/eeglab14_1_2b')
addpath(genpath('/project/def-gnorthof/sorenwt/MATLAB'))
rmpath(genpath('/project/def-gnorthof/sorenwt/MATLAB/osl-core-master'))


settings = SetDefaults(settings);

cd(settings.inputdir)

if ~isempty(find(contains(settings.steps,'tf_filter')))
    disp('Performing time-frequency analysis...')
    tic;
    settings = TTV_ERSP_tf_func(settings);
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
    TTV_ERSP_results_func(settings)
    time_post = toc;
    disp(['Calculation took ' num2str(time_post) ' seconds'])
end

if ~isempty(find(contains(settings.steps,'restcorr')))
    disp('Calculating and correlating resting state parameters...')
    tic;
    TTV_ERSP_rest_func(settings)
    
    time_post = toc;
    disp(['Resting state analysis took ' num2str(time_post) ' seconds'])
end

if ~isempty(find(contains(settings.steps,'behav')))
    disp('Calculating behaviour relationships...')
    tic;
   TTV_ERSP_behaviour_func(settings) 
    time_post = toc;
    disp(['Behaviour took ' num2str(time_post) ' seconds'])
end

if ~isempty(find(contains(settings.steps,'FDR')))
    disp('Correcting across analyses with FDR...')
    tic;
   TTV_ERSP_FDR_correct(settings)
   time_post = toc;
   disp(['FDR correction took ' num2str(time_post) ' seconds'])
end

if ~isempty(find(contains(settings.steps,'figures')))
    disp('Generating figures...')
    tic;
    system(['mkdir ' settings.outputdir '/' settings.datasetname '_figures'])
    TTV_ERSP_figures_func(settings)
    time_post = toc;
    disp(['Creating figures took ' num2str(time_post) ' seconds'])
end

alltime_post = cputime;

disp(['Analysis finished'])
disp(['All steps took ' num2str(alltime_post-alltime_pre) ' seconds'])
disp(['Results saved in ' settings.outputdir])

end


function settingsout = SetDefaults(settingsin)

settingsout = settingsin;

if ~isfield(settingsout,'steps')
   settingsout.steps = {'all'}; 
end

if strcmpi(settingsout.steps{1},'all')
    settingsout.steps = {'tf_filter','results','restcorr','FDR','figures'};
end

if ~isfield(settingsin,'srate')
   settingsout.srate = 200;
end

if ~isfield(settingsin,'units')
   settingsout.units = 'prcchange'; 
end

if ~isfield(settingsin.tfparams,'trials')
    settingsout.tfparams.trials = 'all';
end

if ~isfield(settingsin.tfparams,'method')
    settingsout.tfparams.method = 'hilbert';
end

if ~isfield(settingsin,'fdr')
   settingsout.fdr = 'yes'; 
end

if ~isfield(settingsin,'load_allmeas')
settingsout.load_allmeas = 'yes';
end

% if isfield(settingsin,'rest') && ~isfield(settingsin.rest,'corrmethod')
%    settingsout.rest.corrmethod = 'subject'; 
% end

settingsout.nfreqs = length(settingsout.tfparams.fbandnames);
settingsout.nbchan = length(settingsout.datasetinfo.label);
end
