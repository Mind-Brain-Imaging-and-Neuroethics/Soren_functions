function [EEG,ftdata] = EEG_fromCSV(filein,columnLabels)

if isstr(columnLabels)
   columnLabels = readtable(columnLabels); 
end

data = readtable(filein);
data.Properties.VariableNames = columnLabels.Properties.VariableNames;

ftdata = struct; 
ftdata.trial{1} = data{:,5:end}';
ftdata.fsample = 1024;

cfg = []; cfg.event = 1:3072:(height(data)-3071); cfg.epoch = [0 3072]; cfg.unit = 'samples';
ftdata = ft_epoch(cfg,ftdata);

for c = 1:length(ftdata.trial)
   ftdata.time{c} = linspace(-1.5,1.5,3072);
   ftdata.trialinfo(c) = data.condition(cfg.event(c));
end
trlinfo = ftdata.trialinfo;

EEG = ft2eeglab(ftdata);

for c = 1:70
    EEG.chanlocs(c).labels = data.Properties.VariableNames{4+c};
end

eegdir = extractBefore(which('eeglab'),'eeglab.m');

EEG = pop_chanedit(EEG,'lookup',fullfile(char(eegdir),'plugins','dipfit2.3','standard_BESA','standard-10-5-cap385.elp'),'eval','chans = pop_chancenter( chans, [],[]);');

EEG = pop_select(EEG,'nochannel',[28 65:69]);

EEG = pop_resample(EEG,512);

ftdata = eeglab2fieldtrip(EEG,'preprocessing','none');
ftdata.trialinfo = trlinfo;