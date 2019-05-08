function [stats] = EasyClusterCorrect_mixedeff(data,mxdesign,datasetinfo,formula,parnum)

%data should be a 1x2 cell array containing the independent and dependent
%variables 
%each cell should be a matrix in the form channels x all observations
%datasetinfo should be a fieldtrip structure which gives the channel
%labels and gradiometer/electrode positions
%pass in the design info as a separate argument
%formula should use the variable names 'ind', 'dep', and 'subs'
%parnum is the size of the parallel pool to use (will be too slow otherwise)

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

cfg = []; cfg.method = 'montecarlo'; cfg.statistic = 'ft_statfun_mixedeff';
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
if size(mxdesign,1) == 1
    mxdesign = mxdesign';
end
cfg.mxdesign = mxdesign;
cfg.mxformula = formula;
cfg.parpool = parnum;
%cfg.mediator = mediator;

stats = ft_timelockstatistics(cfg,tlock{1},tlock{2});