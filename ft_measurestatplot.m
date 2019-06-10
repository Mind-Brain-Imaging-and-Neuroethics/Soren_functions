function figs = ft_measurestatplot(cfg,data,stats)
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
%      measname: names of the measures, to plot as figure titles (default =
%         from measure handles)
%      lay: a fieldtrip layout if using meg data
%      plotmode: 'topo','violin', or 'combined'. Plots either the topography 
%         of the measures or a violin plot of the means for each group. 
%         'combined' plots the violin plots below and the topoplots above
%         (default = 'combined')
%      measname: a cell array indicating the names of the measures (default
%         = none)
%      colormap: the colormap to use when plotting topographies (default =
%         parula)
%      datasetinfo: used for source-space plotting. Should contain the
%         fields:
%         
%         atlas: the atlas used to parcellate the data
%         sourcemodel: a head shape or source model that corresponds to the
%         atlas
%      savefig: save the figures as Matlab figures (default = 'no')
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

cfg = setdefault(cfg,'measname',cellfun(@func2str,data{1}.meas(cfg.meas),'UniformOutput',false));

cfg = setdefault(cfg,'colormap',parula);

if ~cfgcheck(cfg,'plotmode')
    cfg.plotmode = 'combined';
end

cfg = setdefault(cfg,'savefig','no');

if nargin < 3
    stats = [];
end

%% Plotting topos

if cfgcheck(cfg,'plotmode','topo')
    if ~isempty(stats)
        if cfgcheck(cfg,'datatype','eeg')
            for c = cfg.meas
                figs(c) = figure;
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    topoplot(mean(data{cc}.data(:,:,c),1),data{1}.chanlocs);
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                    colormap(cfg.colormap)
                end
                
                subplot(1,length(data)+1,length(data)+1)
                if isfield(stats{c},'cluster')
                    cluster_topoplot(mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chanlocs,stats{c}.p,stats{c}.cluster.prob < 0.05)
                elseif isfield(stats{c},'fdr')
                    cluster_topoplot(mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chanlocs,stats{c}.p,stats{c}.fdr < 0.05)
                else
                    cluster_topoplot(mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chanlocs,stats{c}.p,stats{c}.p < 0.05)
                end
                title([cfg.cond{1} ' - ' cfg.cond{2}])
                FixAxes(gca,16)
                colormap(cfg.colormap)
                Normalize_Clim(gcf)
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = cfg.measname{c}; cbar.Label.FontSize = 14;
                set(figs(c),'name',cfg.measname{c})
            end
        elseif cfgcheck(cfg,'datatype','meg')
            for c = cfg.meas
                figs(c) = figure;
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    ft_topoplot_vec(cfg.lay,mean(data{cc}.data(:,:,c),1),data{1}.chan);
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                    colormap(cfg.colormap)
                end
                
                subplot(1,length(data)+1,length(data)+1)
                if isfield(stats{c},'cluster')
                    ft_cluster_topoplot(cfg.lay,mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chan,stats{c}.p,stats{c}.cluster.prob < 0.05)
                elseif isfield(stats{c},'fdr')
                    ft_cluster_topoplot(cfg.lay,mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chan,stats{c}.p,stats{c}.fdr < 0.05)
                else
                    ft_cluster_topoplot(cfg.lay,mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        data{1}.chan,stats{c}.p,stats{c}.p < 0.05)
                end
                title([cfg.cond{1} ' - ' cfg.cond{2}])
                Normalize_Clim(gcf)
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = cfg.measname{c}; cbar.Label.FontSize = 14;
                FixAxes(gca,16)
                colormap(cfg.colormap)
                                set(figs(c),'name',cfg.measname{c})

            end
            
        elseif cfgcheck(cfg,'datatype','source')
            for c = cfg.meas
                figs(c) = figure;
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    ft_cluster_sourceplot(mean(data{cc}.data(:,:,c),1),cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,...
                        ones(size(mean(data{cc}.data(:,:,c),1))));
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                    colormap(cfg.colormap)
                end
                
                subplot(1,length(data)+1,length(data)+1)
                if isfield(stats{c},'cluster')
                    ft_cluster_sourceplot(mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,stats{c}.cluster.prob < 0.05)
                elseif isfield(stats{c},'fdr')
                    ft_cluster_topoplot(mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,stats{c}.fdr < 0.05)
                else
                    ft_cluster_topoplot(cfg.lay,mean(data{1}.data(:,:,c),1)-mean(data{2}.data(:,:,c),1),...
                        cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,stats{c}.p < 0.05)
                end
                title([cfg.cond{1} ' - ' cfg.cond{2}])
                Normalize_Clim(gcf)
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = cfg.measname{c}; cbar.Label.FontSize = 14;
                %FixAxes(gca,16)
                colormap(cfg.colormap)
                                set(figs(c),'name',cfg.measname{c})

            end
            
        end
    else
        if cfgcheck(cfg,'datatype','eeg')
            for c = cfg.meas
                figs(c) = figure;
                for cc = 1:length(data)
                    subplot(1,length(data),cc)
                    topoplot(mean(data{cc}.data(:,:,c),1),data{1}.chanlocs);
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                    colormap(cfg.colormap)
                end
                Normalize_Clim(gcf)
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = cfg.measname{c}; cbar.Label.FontSize = 14;
                            set(figs(c),'name',cfg.measname{c})

            end
            
        elseif cfgcheck(cfg,'datatype','meg')
            for c = cfg.meas
                figs(c) = figure;
                for cc = 1:length(data)
                    subplot(1,length(data),cc)
                    ft_topoplot_vec(cfg.lay,mean(data{cc}.data(:,:,c),1),data{1}.chan);
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                    colormap(cfg.colormap)
                end
                Normalize_Clim(gcf)
                cbar = colorbar('Location','eastoutside'); cbar.Label.String = cfg.measname{c}; cbar.Label.FontSize = 14;
                            set(figs(c),'name',cfg.measname{c})

            end
        elseif cfgcheck(cfg,'datatype','source')
            for c = cfg.meas
                figs(c) = figure;
                for cc = 1:length(data)
                    subplot(1,length(data)+1,cc)
                    ft_cluster_sourceplot(mean(data{cc}.data(:,:,c),1),cfg.datasetinfo.sourcemodel,cfg.datasetinfo.atlas,...
                        ones(size(mean(data{cc}.data(:,:,c),1))));
                    title(cfg.cond{cc})
                    FixAxes(gca,16)
                    colormap(cfg.colormap)
                end
                            set(figs(c),'name',cfg.measname{c})
    
            end
        end
    end
elseif cfgcheck(cfg,'plotmode','violin')
    %% Plotting violins
    for i = cfg.meas
        figs(i) = figure;
        for c = 1:length(data)
            datastruct.(cfg.cond{c}) = mean(data{c}.data(:,:,i),2);
        end
        ylabel(cfg.measname{i})
        violinplot(datastruct)
        FixAxes(gca,16)
        set(figs(i),'name',cfg.measname{c});
    end
elseif cfgcheck(cfg,'plotmode','combined')
    for i = cfg.meas
        tmpcfg = cfg; tmpcfg.plotmode = 'topo'; tmpcfg.meas = i;
        topofig = ft_measurestatplot(tmpcfg,data,stats);
        topofig = topofig(i);
        
        tmpcfg.plotmode = 'violin'; tmpcfg.meas = i;
        violinfig = ft_measurestatplot(tmpcfg,data,stats);
        violinfig = violinfig(i);
        
        figs(i) = figure;
        p = panel('no-manage-font');
        p.pack('v',{40 60})
        nplots = length(data)+(~isempty(stats));
        p(1).pack('h',[repmat({96/nplots},1,nplots) {4}]);
        figaxes = findobj('Parent',topofig,'Type','axes');
        for c = 1:nplots
            p(1,c).select(figaxes(nplots-c+1));
        end
        %cbar = findobj('Parent',topofig,'Type','colorbar');
        %p(1,length(data)+1).select(cbar);
        colormap(cfg.colormap)
        
        
        figaxes = findobj('Parent',violinfig,'Type','axes');
        p(2).select(figaxes)
        
        p.margintop = 10;
        p.marginleft = 20;
        p(1).marginbottom = 5;
        
        AddFigureLabel(p(1,1).axis,'A','yes')
        AddFigureLabel(p(2).axis,'B')
        
        set(figs(i),'name',cfg.measname{c})
        
        close(topofig)
        close(violinfig)
    end
end

if cfgcheck(cfg,'savefig','yes')
    for c = 1:length(figs)
       savefig(figs(c),['Fig ' num2str(c) '-' cfg.measname{c} '.fig']); 
    end
end