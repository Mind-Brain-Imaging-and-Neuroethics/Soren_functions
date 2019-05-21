function ft_measurestatplot(cfg,data,stats)
% ft_measurestatplot plots the output of ft_measurestatistics or
% ft_applymeasure (if stats input is empty)
%
% Inputs:
%
% cfg: a structure with the following fields:
%      datatype: 'eeg', 'meg', or 'source' (default = 'eeg')
%      cond: names for the conditions to be plotted. Should be a cell
%         array of the same length as data (default = {'Condition
%         1','Condition 2', etc}
%      channel: a channel input that works with ft_channelselection
%         (default = 'all')
%      meas: indices of the measures you want to plot (default = 'all' -
%         makes a new figure for each measure
%      lay: a fieldtrip layout if using meg data
%      plotmode: 'topo' or 'violin'. Plots either the topography of the
%         measures or a violin plot of the means for each group (default =
%         'topo')
%
% data: a cell array of outputs structs from ft_applymeasure
%
% stats: a stats cell array from ft_measurestatistics. If no input is
%      given, only the topography of the data is plotted

%% Set up defaults
if ~cfgcheck(cfg,'datatype')
    cfg.datatype = 'eeg';
end

if ~cfgcheck(cfg,'cond')
    cfg.cond = repmat({'Condition'},1,length(data));
    for c = 1:length(cfg.cond)
        cfg.cond{c} = [cfg.cond{c} '_' num2str(c)];
    end
end

if ~cfgcheck(cfg,'channel')
    cfg.channel = 'all';
end

if ~cfgcheck(cfg,'meas')
    cfg.meas = 1:length(data{1}.meas);
end

if ~cfgcheck(cfg,'plotmode')
   cfg.plotmode = 'topo'; 
end

%% Plotting topos

if cfgcheck(cfg,'plotmode','topo')
    if exist('stats','var')
        if cfgcheck(cfg,'datatype','eeg')
            for c = 1:length(cfg.meas)
                figure
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    topoplot(mean(data{cc}.data(:,:,c),1),data{1}.chanlocs);
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                end
                
                subplot(1,length(data)+1,length(data)+1)
                if isfield(stats{c},'cluster')
                    cluster_topoplot(mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chanlocs,stats{c}.p,stats{c}.cluster.mask)
                elseif isfield(stats{c},'fdr')
                    cluster_topoplot(mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chanlocs,stats{c}.p,stats{c}.fdr < 0.05)
                else
                    cluster_topoplot(mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chanlocs,stats{c}.p,stats{c}.p < 0.05)
                end
                title([cfg.cond{1} ' - ' cfg.cond{2}])
                FixAxes(gca,16)
                Normalize_Clim(gcf)
                colorbar('EastOutside')
                
            end
        elseif cfgcheck(cfg,'datatype','meg')
            for c = 1:length(cfg.meas)
                figure
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    ft_topoplot_vec(cfg.lay,mean(data{cc}.data(:,:,c),1),data{1}.chan);
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                end
                
                subplot(1,length(data)+1,length(data)+1)
                if isfield(stats{c},'cluster')
                    ft_cluster_topoplot(cfg.lay,mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chan,stats{c}.p,stats{c}.cluster.mask)
                elseif isfield(stats{c},'fdr')
                    ft_cluster_topoplot(cfg.lay,mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chan,stats{c}.p,stats{c}.fdr < 0.05)
                else
                    ft_cluster_topoplot(cfg.lay,mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chan,stats{c}.p,stats{c}.p < 0.05)
                end
                title([cfg.cond{1} ' - ' cfg.cond{2}])
                Normalize_Clim(gcf)
                colorbar('EastOutside')
                FixAxes(gca,16)
            end
            
        elseif cfgcheck(cfg,'datatype','source')
            
        end
    else
        if cfgcheck(cfg,'datatype','eeg')
            for c = 1:length(cfg.meas)
                figure
                for cc = 1:length(data)
                    subplot(1,length(data),cc)
                    topoplot(mean(data{cc}.data(:,:,c),1),data{1}.chanlocs);
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                end
                Normalize_Clim(gcf)
                colorbar('EastOutside')
            end
            
        elseif cfgcheck(cfg,'datatype','meg')
            for c = 1:length(cfg.meas)
                figure
                for cc = 1:length(data)
                    subplot(1,length(data),cc)
                    ft_topoplot_vec(cfg.lay,mean(data{cc}.data(:,:,c),1),data{1}.chan);
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                end
                Normalize_Clim(gcf)
                colorbar('EastOutside')
            end
            
            
        end
    end
elseif cfgcheck(cfg,'plotmode','violin')
    %% Plotting violins
    for i = 1:length(cfg.meas)
        figure
        for c = 1:length(data)
            datastruct.(cfg.cond{c}) = mean(data{c}.data(:,:,i),2);
        end
        violinplot(datastruct)
    end
end