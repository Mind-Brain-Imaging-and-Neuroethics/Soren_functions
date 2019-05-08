function [cont_data] = RAM_preproc_continuous(subid,cont_data)

%% Filter before epoching because of the low frequencies
% Also filter first so line noise is dealt with before bad channels
cfg = []; cfg.lpfilter = 'yes'; cfg.lpfreq = [200];
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
D = spm_eeg_ft2spm(cont_data,fullfile(pwd,subid,'tmp3'));
D.elec = cont_data.elec;
S.dataset = fullfile(pwd,subid,'tmp3');
S.FileOut = fullfile(pwd,subid,'tmp4');  % If ommited, the input file will be updated, no copy of the file will be created
S.trainBase = '/home/soren/Documents/MATLAB/ImaGIN2/toolbox';   % Directory where ImaGIN_trainBaseFeatures.mat is located
[D,bIdx] = ImaGIN_BadChannel(S,D);
%!rm tmp.mat tmp.dat tmp2.mat tmp2.dat tmp2_bIdx.txt tmp2_bChans.txt recordings_monopolar_tmp2.txt
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

%% Preprocess the continuous data 

% Epoch into arbitrary 2s segments

cfg = []; cfg.event = 1:2*cont_data.fsample:floor(length(cont_data.time{1})); 
cfg.event(end) = []; cfg.epoch = [0 2*cont_data.fsample-1];
cont_data = ft_epoch(cfg,cont_data);

% Artifact rejection

cont_data = ft_autoartifact(cont_data,{'jump'},0);

cont_data = ft_struct2single(cont_data);











