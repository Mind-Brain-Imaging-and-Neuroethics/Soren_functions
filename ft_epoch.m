function [trialdata] = ft_epoch(cfg,contdata)

%inputs in cfg:
% event: times of event onset
% epoch: size of epoch around the event
% unit: 'samples' or 'seconds'

if ~isfield(cfg,'unit')
   cfg.unit = 'samples'; 
end

trialdata = contdata;
trialdata.time = [];
trialdata.trial = [];

if strcmpi(cfg.unit,'seconds')
   cfg.event = FindClosest(contdata.time,cfg.event);
   cfg.epoch = cfg.epoch*contdata.fsample;
end


for c = 1:length(cfg.event)
     trialdata.trial{c} = contdata.trial{1}(:,...
         (cfg.event(c)+cfg.epoch(1)):(cfg.event(c)+cfg.epoch(2)-1));
    trialdata.time{c} = linspace(cfg.epoch(1)/contdata.fsample,cfg.epoch(2)/contdata.fsample,size(trialdata.trial{c},2));
    trialdata.sampleinfo(c,:) = [cfg.event(c)+cfg.epoch(1) cfg.event(c)+cfg.epoch(2)-1];
end