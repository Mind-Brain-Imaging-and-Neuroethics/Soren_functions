function [data,cont_data] = RAM_preproc(subid,cont_data,latencies,continuous)

%% Filter before epoching because of the low frequencies
% Also filter first so line noise is dealt with before bad channels

cfg = []; cfg.bpfilter = 'yes'; cfg.bpfreq = [0.03 200]; cfg.bpinstabilityfix = 'split';
cfg.bsfilter = 'yes'; cfg.bsfreq = [57 63; 118 122; 178 182];
cont_data = ft_preprocessing(cfg,cont_data);


%% Find and remove bad channels
% Reference: Tuyisenge et al. 2018. Automatic bad channel detection in intracranial 
% electroencephalographic recordings using ensemble machine learning. Clinical Neurophysiology

for c = 1:length(cont_data.elec.label)
    cont_data.elec.chantype(c) = {'EEG'};
end
cont_data.elec.type = 'EEG';
cont_data.hdr.chantype = cont_data.elec.chantype;
cont_data.hdr.label = cont_data.label;

spm
D = spm_eeg_ft2spm(cont_data,fullfile(pwd,subid,'tmp'));
D.elec = cont_data.elec;
S.dataset = fullfile(pwd,subid,'tmp');
S.FileOut = fullfile(pwd,subid,'tmp2');  % If ommited, the input file will be updated, no copy of the file will be created
S.trainBase = '/home/soren/Documents/MATLAB/ImaGIN2/toolbox';   % Directory where ImaGIN_trainBaseFeatures.mat is located
[D,bIdx] = ImaGIN_BadChannel(S,D);
!rm tmp.mat tmp.dat tmp2.mat tmp2.dat tmp2_bIdx.txt tmp2_bChans.txt recordings_monopolar_tmp2.txt
!rm -r ScreenShot

D = [];

ft_defaults

cfg = []; cfg.channel = cont_data.label; cfg.channel(bIdx) = [];
cont_data = ft_selectdata(cfg,cont_data);

%% Rereference to average

cfg             = [];
cfg.channel     = {'all'};
cfg.reref       = 'yes';
cfg.refchannel  = 'all';
cfg.refmethod   = 'avg';
cont_data = ft_preprocessing(cfg, cont_data);

%% Epoch data
cfg = []; cfg.event = latencies; cfg.epoch = [-3*cont_data.fsample 3*cont_data.fsample];
data = ft_epoch(cfg,cont_data);

% Code for bipolar montage if that is needed later
%
% depths = unique(erase(data.label,cellstr(num2str([0:9]'))));
% for d = 1:numel(depths)
%     currdepths = [depths{d} '*'];
%     cfg            = [];
%     cfg.channel    = ft_channelselection(currdepths, data.label);
%     out = [];
%     for c = 1:length(cfg.channel)
%        if isempty(str2num(cfg.channel{c}(length(depths{d})+1)))
%           out = [out c]; 
%        end
%     end
%     cfg.channel(out) = [];
%     cfg.reref      = 'yes';
%     cfg.refchannel = 'all';
%     cfg.refmethod  = 'bipolar';
%     cfg.updatesens = 'yes';
%     reref_depths{d} = ft_preprocessing(cfg, data);
% end
% 
% % Concatenate data
% cfg            = [];
% cfg.appendsens = 'yes';
% data = ft_appenddata(cfg, reref_depths{:});

%% Artifact rejection
% Reference: Defaults from Fieldtrip examples

data = ft_autoartifact(data,{'jump'},0);

data = ft_struct2single(data);

%% Preprocess the continuous data 
if continuous

% Epoch into arbitrary 2s segments

cfg = []; cfg.event = 1:2*cont_data.fsample:floor(length(cont_data.time{1})); 
cfg.event(end) = []; cfg.epoch = [0 2*cont_data.fsample-1];
cont_data = ft_epoch(cfg,cont_data);

% Artifact rejection

cont_data = ft_autoartifact(cont_data,{'jump'},0);

cont_data = ft_struct2single(cont_data);
end











