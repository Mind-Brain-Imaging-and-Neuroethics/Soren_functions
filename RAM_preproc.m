function [data] = RAM_preproc(subid,cont_data,latencies,continuous)

%% Filter before epoching because of the low frequencies
% Also filter first so line noise is dealt with before bad channels

cfg = [];
if continuous
    cfg.lpfilter = 'yes'; cfg.lpfreq = 200;
else
    cfg.bpfilter = 'yes'; cfg.bpfreq = [0.03 200]; cfg.bpinstabilityfix = 'split';
end
cfg.bsfilter = 'yes'; cfg.bsfreq = [57 63; 118 122; 178 182];
cont_data = ft_preprocessing(cfg,cont_data);

%% Find and remove bad channels
% Reference: Tuyisenge et al. (2018). Automatic bad channel detection in intracranial
% electroencephalographic recordings using ensemble machine learning. Clinical Neurophysiology.

for c = 1:length(cont_data.elec.label)
    cont_data.elec.chantype(c) = {'EEG'};
end
cont_data.elec.type = 'EEG';
cont_data.hdr.chantype = cont_data.elec.chantype;
cont_data.hdr.label = cont_data.label;

spm
addpath([toolboxdir('signal') '/signal'])
D = spm_eeg_ft2spm(cont_data,fullfile(pwd,subid,'tmp'));
D.elec = cont_data.elec;
S.dataset = fullfile(pwd,subid,'tmp');
S.FileOut = fullfile(pwd,subid,'tmp2');  % If ommited, the input file will be updated, no copy of the file will be created
S.trainBase = '/home/soren/Documents/MATLAB/ImaGIN2/toolbox';   % Directory where ImaGIN_trainBaseFeatures.mat is located
[D,bIdx] = ImaGIN_BadChannel(S,D);
system(['rm ' fullfile(pwd,subid,'tmp.mat') ' ' fullfile(pwd,subid,'tmp.dat') ' ' fullfile(pwd,subid,'tmp2.mat')...
    fullfile(pwd,subid,'tmp2.dat') ' ' fullfile(pwd,subid,'tmp2_bIdx.txt') ' ' ...
    fullfile(pwd,subid,'tmp2_bChans.txt') ' ' fullfile(pwd,subid,'recordings_monopolar_tmp2.txt')])
system('rm -r ScreenShot')

D = [];

ft_defaults

cfg = []; cfg.channel = cont_data.label; cfg.channel(bIdx) = [];
cont_data = ft_selectdata(cfg,cont_data);

%% Epoch data
if continuous
    % Epoch into arbitrary 2s segments
    cfg = []; cfg.event = 1:2*cont_data.fsample:floor(length(cont_data.time{1}));
    cfg.event(end) = []; cfg.epoch = [0 2*cont_data.fsample-1];
    data = ft_epoch(cfg,cont_data);
    cont_data = [];
else
    % Epoch according to provided latencies
    cfg = []; cfg.event = latencies; cfg.epoch = [-1.5*cont_data.fsample 1.5*cont_data.fsample];
    data = ft_epoch(cfg,cont_data);
    cont_data = [];
end

data.trialinfo = ones(length(data.sampleinfo),1);

%% Resample to 500 Hz

cfg = []; cfg.resamplefs = 500; 
cont_data = ft_resampledata(cfg,cont_data);

%% Convert to MNE, get rejection threshold with autoreject
% Reference: Jas et al. (2016). Autoreject: Automated artifact rejection
% for MEG and EEG data. Neuroimage.

save(fullfile(pwd,subid,[subid '_epochs.mat']),'data')
pyscript = fopen([subid '_pyscript.py'],'w')
fprintf(pyscript,'import sys \n')
fprintf(pyscript,'sys.path.insert(0, ''/home/soren/Documents/MATLAB/Functions'') \n')
fprintf(pyscript,'from mne_preproc import autoreject_threshold \n')
fin = fullfile(pwd,subid,[subid '_epochs.mat']);
fout = fullfile(pwd,subid,[subid '_threshold.json']);
fprintf(pyscript,['autoreject_threshold(''' fin ...
    ''',''' fout ''')'])
[a,b] = system(['python ' subid '_pyscript.py']);
system(['rm ' subid '_pyscript.py'])

if a == 1
   error(['Error using autoreject - message was: ' b])
end

threshold = jsonread(fullfile(pwd,subid,[subid '_threshold.json']));
threshold = threshold.eeg;

%% Reject trials over threshold, find more bad channels

badtrl = zeros(length(data.label),length(data.trial));
badchans = zeros(length(data.label),1);
for c = 1:length(data.label)
    for cc = 1:length(data.trial)
        pk2pk = max(data.trial{cc}(c,:)) - min(data.trial{cc}(c,:));
        badtrl(c,cc) = pk2pk > threshold;
    end
end

if continuous
    badprct = 0.1;
else
    badprct = 0.2;
end

% for c = 1:length(data.label)
%     if sum(badtrl(c,:)) > (length(data.trial)*badprct)
%         badchans(c) = 1;
%         badtrl(c,:) = zeros(1,length(data.trial));
%     end
% end

rejtrl = sum(badtrl,1) > 0;
while sum(rejtrl) > badprct*length(data.trial)
    [~,c] = max(sum(badtrl,2));
    badchans(c) = 1;
    badtrl(c,:) = zeros(1,length(data.trial));
    rejtrl = sum(badtrl,1) > 0;
end

badtrl = sum(badtrl,1) > 0;

% if sum(badtrl) > (length(badtrl)/2)
%     data = NaN;
%     error('Too many bad trials - excluding subject')
% end

if sum(badchans) > 0.5*length(data.label)
    data = NaN;
    error('Too many bad trials - excluding subject')
end

cfg = []; cfg.trials = find(~badtrl); cfg.channel = find(~badchans);
data = ft_selectdata(cfg,data);

%% Rereference to average

cfg             = [];
cfg.channel     = {'all'};
cfg.reref       = 'yes';
cfg.refchannel  = 'all';
cfg.refmethod   = 'avg';
data = ft_preprocessing(cfg, data);


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











