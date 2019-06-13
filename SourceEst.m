function [roidata,voxeldata,sources] = SourceEst_MEG(data,headmodel,sourcemodel,atlas,noisecov,opts)
% opts should have the following fields
%    atlasparam: the field of the atlas containing the parcellation
%      information (default = 'parcellation')
%    atlaslabel: the field of the atlas containing the labels for each
%      region (default = 'parcellationlabel')
%    interp: interpolate the sourcemodel onto the atlas ('yes' or 'no' -
%      default = 'no')
%    datatype: 'EEG' or 'MEG' - for channel selection (default = 'MEG')


%[~,ftpath] = ft_version;
%ftpath = '/group/northoff/share/fieldtrip-master';

% % Convert the data to Fieldtrip format
% data = eeglab2fieldtrip(EEG,'preprocessing','none');

if nargin < 6
    opts = struct;
end

opts = setdefault(opts,'atlasparam','parcellation');
opts = setdefault(opts,'atlaslabel','parcellationlabel');
opts = setdefault(opts,'interp','no');
opts = setdefault(opts,'datatype','MEG');

% Load headmodel, sourcemodel, atlas

if isstr(headmodel)
    load(headmodel) %variable is "vol"
end

if isstr(sourcemodel)
    sourcemodel = ft_read_headshape(sourcemodel);
end

if isstr(atlas)
    atlas = ft_read_atlas(atlas);
end

% % Interpolate template surface on atlas
if strcmpi(opts.interp,'yes')
    cfg = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter = opts.atlasparam;
    sourcemodel = ft_sourceinterpolate(cfg, atlas, sourcemodel);
    atlas.(opts.atlasparam) = sourcemodel.(opts.atlasparam);
end

% Create forward model
cfg = []; 
if strcmpi(opts.datatype,'MEG')
cfg.grad = data.grad; 
elseif strcmpi(opts.datatype,'EEG')
cfg.elec = data.elec
end
cfg.channel = {opts.datatype};
sourcemodel = ft_convert_units(sourcemodel,'mm');
headmodel = ft_convert_units(headmodel,'mm');
cfg.grid.pos = sourcemodel.pos; cfg.grid.inside = 1:size(sourcemodel.pos,1);
cfg.headmodel = headmodel;
leadfield = ft_prepare_leadfield(cfg);

% % Select only the relevant channels
cfg = []; cfg.channel = {opts.datatype};
data = ft_selectdata(cfg,data);

% If continuous, epoch into arbitrary 2-second segments
if length(data.trial) == 1
    cfg = []; cfg.event = 1:2*data.fsample:length(data.time{1}); cfg.epoch = [0 (2*data.fsample)-1];
    cfg.event(end) = [];
    data = ft_epoch(cfg,data);
end

cfg = []; cfg.keeptrials = 'yes';
tlock = ft_timelockanalysis(cfg,data);
tlock.cov = noisecov;

% Estimate sources
cfg = [];
cfg.method = 'eloreta';
cfg.grid = leadfield;
cfg.eloreta.keepfilter = 'yes';
cfg.eloreta.normalize = 'yes';
cfg.headmodel = headmodel;
sources = ft_sourceanalysis(cfg,tlock);

% Get voxel time courses
voxeldata = struct;
datacat = cat(2,data.trial{:});
for c = 1:size(sources.pos,1)
    tmp = sources.avg.filter{c}*datacat;
    [u,s,v] = svd(tmp,'econ');
    tmp = u(:,1)'*sources.avg.filter{c}*datacat;
    source_datacat(c,:) = tmp;
end
voxeldata.trial{1} = source_datacat;
voxeldata.time{1} = linspace(0,length(voxeldata.trial{1}/data.fsample),1/data.fsample);
voxeldata.label = cellstr(num2str([1:size(sources.pos,1)]'));
voxeldata.fsample = data.fsample;
source_datacat = []; %saving memory

% Get ROI time courses
roidata = voxeldata;
roidata.label = atlas.(opts.atlaslabel);
roidata.trial = [];
roidata.trial{1} = NaN(length(roidata.label),length(voxeldata.trial{1}));
for cc = 1:length(roidata.label)
    roidata.trial{1}(cc,:) = mean(voxeldata.trial{1}(find(atlas.(opts.atlasparam) == cc),:),1);
end
%roidata.trial{1}(find(isnan(roidata.trial{1}(:,1))),:) = [];

end
