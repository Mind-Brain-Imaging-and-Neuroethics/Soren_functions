function [freqdata,specs] = IRASA_tf(cfg,data,specs)
% cfg:
%     oscifrac: output oscillatory or fractal spectrum
%     winsize: window size in seconds
%     modifywindow: 'yes' or 'no' - if extentsion changes window size by 
%         less than 15%, extends window to optimize FFTs (default = 'yes')
%     toi: times of interest to center the window on
%     foi: optional cell array defining frequency bands of interest. Each
%         cell should be a two-element vector of the form [min_freq
%         max_freq] (default: outputs all frequencies from IRASA)
%     hset: hset for irasa resampling (default = [1.1:0.05:1.95
%     2.05:0.05:2.9])
%     parflag: use parallel pool - 'yes' or 'no' 
% data: a standard fieldtrip data structure (output from ft_preprocessing)
% specs: optional argument - allows inputing the raw spectra from a
%    previous iteration in order to get the time-frequency data in a
%    fieldtrip-like format

% Outputs:
%
% tfout: time-frequency decomposition in the form trials x chans x time x
%    freqs
% specsout: the raw spectra from the IRASA decomposition

addpath([toolboxdir('signal') '/signal'])

if ~isfield(cfg,'modifywindow')
    cfg.modifywindow = 'yes';
end

if ~isfield(cfg,'hset')
    cfg.hset =  [1.1:0.05:1.95 2.05:0.05:2.9];
end

cfg = setdefault(cfg,'parflag','no');

if strcmpi(cfg.modifywindow,'yes')
    winlensamples = cfg.winsize*data.fsample;
    sampspow2 = (2^nextpow2(winlensamples))/0.9;
    if (sampspow2-winlensamples)/winlensamples < 0.15
        cfg.winsize = (sampspow2+1)/data.fsample;
        disp(['New window size is ' num2str(cfg.winsize)])
    end
end

if nargin < 3
   calcspecs = 1;
else
    calcspecs = 0;
end

if strcmpi(cfg.parflag,'yes')
    freqdata = cell(1,length(data.trial));
else
    freqdata = struct;
end

if calcspecs
    specs = cell(1,length(data.trial));
end
trials = data.trial;
times = data.time;
fsample = data.fsample;

if strcmpi(cfg.parflag,'yes')
    parfor i = 1:length(data.trial)
        datawindows = getWindows(times{i},cfg.winsize,cfg.toi,trials{i});
        for c = 1:length(datawindows)
	   
            if calcspecs
                specs{i}(c) = amri_sig_fractal(datawindows{c}',fsample,'hset',cfg.hset);
            end
            
            if ~isfield(cfg,'foi')
                freqdata{i}(1,:,:,c) = specs{i}(c).(cfg.oscifrac);
            else
                for cc = 1:length(cfg.foi)
                    findex = intersect(find(specs{i}(c).freq > cfg.foi{cc}(1)),find(specs{i}(c).freq < cfg.foi{cc}(2)));
                    freqdata{i}(1,:,cc,c) = trapz(specs{i}(c).freq(findex),specs{i}(c).(cfg.oscifrac)(findex,:),1);
                end
            end
        end
    end
    tmp = cat(1,freqdata{:});
    freqdata = struct;
    freqdata.fourierspctrm = tmp;
    clear tmp
else
    for i = 1:length(data.trial)
        datawindows = getWindows(data.time{i},cfg.winsize,cfg.toi,data.trial{i});
        for c = 1:length(datawindows)
            if calcspecs
                specs{i}(c) = amri_sig_fractal(datawindows{c}',data.fsample,'hset',cfg.hset);
            end
            
            if ~isfield(cfg,'foi')
                freqdata.fourierspctrm(i,:,:,c) = specs{i}(c).(cfg.oscifrac);
            else
                for cc = 1:length(cfg.foi)
                    findex = intersect(find(specs{i}(c).freq > cfg.foi{cc}(1)),find(specs{i}(c).freq < cfg.foi{cc}(2)));
                    freqdata.fourierspctrm(i,:,cc,c) = mean(specs{i}(c).(cfg.oscifrac)(findex,:),1);
                end
            end
        end
    end
end

if ~isfield(cfg,'foi')
    freqdata.freq = specs{1}(1).freq;
else
    freqdata.freq = cell2mat(cellfun(@mean,cfg.foi,'UniformOutput',false));
end

freqdata.dimord = 'rpttap_chan_freq_time';
freqdata.label = data.label;
freqdata.time = cfg.toi;

if isfield(data,'elec')
    freqdata.elec = data.elec;
elseif isfield(data,'grad')
    freqdata.grad = data.grad;
end
