function [stats] = EasyClusterCorrect_mediation(data,mediator,datasetinfo)

%data should be a 1x2 cell array containing the independent and dependent
%variables - the mediator should be specified as a separate array
%each cell should be a matrix in the form channels x observations
%datasetinfo should be a fieldtrip structure which gives the channel
%labels and gradiometer/electrode positions

% if nargin > 2
%    oldcfg = cfg;
%    cfg = [];
% end

fields = fieldnames(datasetinfo);
for c = 1:2
   tlock{c} = struct;
   tlock{c}.avg = [mean(data{c},2) mean(data{c},2) mean(data{c},2)];
   tlock{c}.trial = cat(3,data{c}',data{c}',data{c}');
   tlock{c}.time = [-1 0 1];
   tlock{c}.dimord = 'rpt_chan_time';
      for q = 1:length(fields)
      tlock{c}.(fields{q}) = datasetinfo.(fields{q});
   end
   %tlock{c}.label = datasetinfo.label;
   %tlock{c}.grad = datasetinfo.grad;
      tlock{c}.fsample = 1;
   tlock{c}.sampleinfo = [[1:size(tlock{c}.trial,1)]' [1:size(tlock{c}.trial,1)]'+3];
end

cfg = []; cfg.method = 'distance'; 
neighbours = ft_prepare_neighbours(cfg,datasetinfo);

cfg = []; cfg.method = 'montecarlo'; cfg.statistic = 'ft_statfun_mediation';
cfg.correctm = 'cluster'; cfg.clusteralpha = 0.05; cfg.clusterstatistic = 'maxsum';
cfg.neighbours = neighbours; cfg.tail = 1; cfg.clustertail = 1; cfg.alpha = 0.05;
cfg.numrandomization = 2000; cfg.spmversion = 'spm12';

design = zeros(1,size(data{1},2) + size(data{2},2));
design(1,1:size(data{1},2)) = 1;
design(1,(size(data{1},2)+1):(size(data{1},2) + size(data{2},2)))= 2;

cfg.design = design; cfg.ivar = 1;

% if isfield(datasetinfo,'grad')
%     cfg.channel = {'MEG'};
% else
%     cfg.channel = {'EEG'};
% end
cfg.channel = {'all'};

cfg.latency = [0 0.5];
cfg.parameter = 'trial';

cfg.mediator = mediator;

stats = ft_timelockstatistics(cfg,tlock{1},tlock{2});