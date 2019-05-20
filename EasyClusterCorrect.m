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
%            'atlasname' (currently only 'aal','mmp', and 'yeo' suppported
%            - for 'yeo' the atlas must be modified to include 'pos' and
%            'tri' fields from an HCP subject's source model)
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
    badlabels = [];

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
    for cc = 1:length(tlock{c}.label)
        if isempty(find(~isnan(tlock{c}.trial(cc,:,1))))
            badlabels = [badlabels tlock{c}.label{cc}];
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
        case 'mmp','yeo'
            
            if strcmpi(datasetinfo.atlasname,'yeo')
                datasetinfo.atlas.parcellationlabel = cellstr(num2str([1:8004]'));
                datasetinfo.atlas.parcellation = datasetinfo.atlas.parcels;
            end
            
            pos = datasetinfo.atlas.pos;
            tri = datasetinfo.atlas.tri;
            
            vox_neighbs = cell(1,length(pos));
            
            for c = 1:length(pos)
                inds = find(tri == c);
                for cc = 1:length(inds)
                    [i1,i2] = ind2sub(size(tri),inds(cc));
                    vox_neighbs{c} = [vox_neighbs{c} atlas.tri(i1,except(1:3,i2))];
                end
                vox_neighbs{c} = unique(vox_neighbs{c});
            end
            
            reg_neighbs = cell(1,length(datasetinfo.atlas.parcellationlabel));
            for c = 1:length(datasetinfo.atlas.parcellationlabel)
                reg_neighbs{c} = cat(2,vox_neighbs{find(datasetinfo.atlas.parcellation == c)});
                for cc = 1:length(reg_neighbs{c})
                    reg_neighbs{c}(cc) = atlas.parcellation(reg_neighbs{c}(cc));
                end
                reg_neighbs{c} = unique(reg_neighbs{c});
                reg_neighbs{c}(find(reg_neighbs{c} == c)) = [];
            end
            
            neighbs = struct;
            for c = 1:length(reg_neighbs)
                neighbs(c).label = datasetinfo.atlas.parcellationlabel{c};
                neighbs(c).neighblabel = {datasetinfo.atlas.parcellationlabel{reg_neighbs{c}}};
            end
            
            vox_neighbs = [];
            reg_neighbs = [];
            pos = [];
            tri = [];
    end
end

for c = 1:2
    cfg = []; cfg.channel = except(1:length(tlock{c}.label),Subset_index(tlock{c}.label,badlabels));
    tlock{c} = ft_selectdata(cfg,tlock{c});
end

rmneighbs = zeros(1,length(neighbs));
for c = 1:length(neighbs)
   if ~isempty(find(strcmpi(neighbs(c).label,badlabels)))
      rmneighbs(c) = 1;
   end
end
neighbs(find(rmneighbs)) = [];

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