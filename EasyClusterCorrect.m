function [stats] = EasyClusterCorrect(data,datasetinfo,statfun,opts)
% EasyClusterCorrect implements fieldtrip's cluster correction algorithm
% across channels only, for an arbitrary measurement
%
% Input arguments: 
%      data: a 1x2 (currently only supports a pairwise contrast) cell array
%         containing the data for testing. Each cell should be a matrix of
%         the form channels x observations
%      datasetinfo: a structure containing the relevant information for
%         clustering. For the following data types, this should be: 
%         eeg: requires fields 'elec' (a fieldtrip electrode structure) and
%            'label'
%         meg: requires fields 'grad' (a fieldtrip gradiometer structure)
%            and 'label'
%         source: currently only region-based correction supported. 
%            Requires fields 'atlas' (a fieldtrip atlas structure) and
%            'atlasname' (currently only 'aal','mmp', and 'yeo' suppported)
%      statfun: the fieldtrip statfun used for the statistics (e.g.
%         ft_statfun_signrank)
%      opts: optional argument which changes some parameters of the
%         clustering. Can have the following fields:
%         nrand: the number of permutations to carry out (default = 2000)
%         minnbchan: the minimum number of neighbours to add a channel to a
%            cluster (default = 1)
%
% Outputs: 
%      stats: a fieldtrip stats structure returned from
%         ft_timelockstatistics

if nargin < 3
    opts = struct;
end

if ~isfield(opts,'nrand')
    opts.nrand = 2000;
end
if ~isfield(opts,'minnbchan')
    opts.minnbchan = 1;
end


fields = fieldnames(datasetinfo);
for c = 1:2
    tlock{c} = struct;
    tlock{c}.avg = [mean(data{c},2) mean(data{c},2) mean(data{c},2)];
    tlock{c}.trial = cat(3,data{c}',data{c}',data{c}');
    tlock{c}.time = [-1 0 1];
    tlock{c}.dimord = 'rpt_chan_time';
    if isfield(datasetinfo,'elec') || isfield(datasetinfo,'grad')
        for q = 1:length(fields)
            tlock{c}.(fields{q}) = datasetinfo.(fields{q});
        end
    elseif isfield(datasetinfo,'atlasname')
        switch datasetinfo.atlasname
            case 'aal'
                tlock{c}.label = datasetinfo.atlas.tissuelabel;
            case 'mmp'
                tlock{c}.label = datasetinfo.atlas.parcellationlabel;
            case 'yeo'
                %do later
        end
    end
    %tlock{c}.label = datasetinfo.label;
    %tlock{c}.grad = datasetinfo.grad;
    tlock{c}.fsample = 1;
    tlock{c}.sampleinfo = [[1:size(tlock{c}.trial,1)]' [1:size(tlock{c}.trial,1)]'+3];
end

if isfield(datasetinfo,'grad') || isfield(datasetinfo,'elec')
    cfg = []; cfg.method = 'distance';
    neighbs = ft_prepare_neighbours(cfg,datasetinfo);
else
    switch datasetinfo.atlasname
        case 'aal'
            tissue = datasetinfo.atlas.tissue;
            tissue(find(tissue == 0)) = NaN;
            
            neighbs = struct;
            for c = 1:length(datasetinfo.atlas.tissuelabel)
                neighbs(c).label = datasetinfo.atlas.tissuelabel{c};
                neighbs(c).neighblabel = cell(1,1);
            end
            
            for c = 1:3
                edges{c} = diff(tissue,1,c);
                edges{c} = edges{c} > 0;
                boundaryindx = find(edges{c});
                for cc = 1:length(boundaryindx)
                    [i1,i2,i3] = ind2sub(size(edges{c}),boundaryindx(cc));
                    indx = {i1 i2 i3};
                    indx{c} = indx{c}+1;
                    roi1val = tissue(i1,i2,i3);
                    roi2val = tissue(indx{:});
                    neighbs(roi1val).neighblabel = [neighbs(roi1val).neighblabel neighbs(roi2val).label];
                    neighbs(roi2val).neighblabel = [neighbs(roi2val).neighblabel neighbs(roi1val).label];
                end
            end
            
            for c = 1:length(neighbs)
                neighbs(c).neighblabel(1) = [];
                neighbs(c).neighblabel = unique(neighbs(c).neighblabel);
            end
        case 'mmp'
            pos = datasetinfo.atlas.pos;
            
            % do later
        case 'yeo'
            % do later
            
    end
end

cfg = []; cfg.method = 'montecarlo'; cfg.statistic = statfun;
cfg.correctm = 'cluster'; cfg.clusteralpha = 0.025; cfg.clusterstatistic = 'maxsum';
cfg.neighbours = neighbs; cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025;
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