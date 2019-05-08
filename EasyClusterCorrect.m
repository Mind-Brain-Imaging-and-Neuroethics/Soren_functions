function [stats] = EasyClusterCorrect(data,datasetinfo,statfun,opts)

%data should be a 1x2 cell array containing the data for testing
%each cell should be a matrix in the form channels x observations
%datasetinfo should be a fieldtrip structure which gives the channel
%labels and gradiometer positions

if nargin < 2
    opts = struct;
end

if ~isfield(opts,'nrand')
    opts.nrand = 2000;
end
if ~isfield(opts,'minnbchan')
    opts.minnbchan = 2;
end


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

cfg = []; cfg.method = 'montecarlo'; cfg.statistic = statfun;
cfg.correctm = 'cluster'; cfg.clusteralpha = 0.025; cfg.clusterstatistic = 'maxsum';
cfg.neighbours = neighbours; cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025;
cfg.numrandomization = opts.nrand; cfg.spmversion = 'spm12'; cfg.minnbchan = opts.minnbchan;

design = zeros(1,size(data{1},2) + size(data{2},2));
design(1,1:size(data{1},2)) = 1;
design(1,(size(data{1},2)+1):(size(data{1},2) + size(data{2},2)))= 2;

cfg.design = design; cfg.ivar = 1;

%if isfield(datasetinfo,'grad')
%    cfg.channel = {'MEG'};
%else
%    cfg.channel = {'EEG'};
%end
cfg.channel = {'all'};

cfg.latency = [0 0.5];
cfg.parameter = 'trial';



stats = ft_timelockstatistics(cfg,tlock{1},tlock{2});