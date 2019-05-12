function TTV_ERSP_figures_func(settings)

fbands = settings.tfparams.fbandnames;

load([settings.outputdir '/' settings.datasetname '_calc.mat'])
load([settings.outputdir '/' settings.datasetname '_results.mat'])
if isfield(settings,'rest')
    load([settings.outputdir '/' settings.datasetname '_restmeas.mat'])
end

if strcmpi(settings.datatype,'EEG')
    %eeglab
end

cd([settings.outputdir '/' settings.datasetname '_figures'])

prestim_pseudo = settings.pseudo.prestim;
prestim_real = settings.real.prestim;
poststim_pseudo = settings.pseudo.poststim;
poststim_real = settings.real.poststim;

%% Figure 1a. Raw TTV, Amplitude, ITC time courses

figure

allrange = min(settings.pseudo.prestim):max(settings.real.poststim);
onset = min(settings.real.poststim);
for c = 1:settings.nfreqs
    subplot(3,settings.nfreqs,c)
    %t = linspace((1/settings.srate)*(min(allrange)-onset),(1/settings.srate)*(max(allrange)-onset),length(allrange));
    t = linspace(0,(1/settings.srate)*length(poststim_pseudo),length(poststim_pseudo));
    plot(t,mean(mean(allmeas{c}.raw.sd(:,poststim_real,:),3),1),'b','LineWidth',1.5)
    hold on
    plot(t,mean(mean(allmeas{c}.raw.sd(:,poststim_pseudo,:),3),1),'b--','LineWidth',1.5)
    
    xlabel('Time (s)')
    ylabel('Across-trial SD')
    title(fbands{c})
    
    subplot(3,settings.nfreqs,c+settings.nfreqs)
    plot(t,mean(mean(allmeas{c}.raw.ersp(:,poststim_real,:),3),1),'r','LineWidth',1.5)
    hold on
    plot(t,mean(mean(allmeas{c}.raw.ersp(:,poststim_pseudo,:),3),1),'r--','LineWidth',1.5)
    
    xlabel('Time (s)')
    ylabel('Hilbert amplitude')
    
    subplot(3,settings.nfreqs,c+2*settings.nfreqs)
    plot(t,mean(mean(allmeas{c}.raw.itc(:,poststim_real,:),3),1),'g','LineWidth',1.5)
    hold on
    plot(t,mean(mean(allmeas{c}.raw.itc(:,poststim_real,:),3),1),'g--','LineWidth',1.5)
    
    xlabel('Time (s)')
    ylabel('Inter-trial coherence')
    
end
%Aesthetics
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)/2 pos(3)*3 pos(4)*2]) %make figure larger

ax = findall(gcf,'Type','Axes');
for c = 1:length(ax)
    set(ax(c),'FontSize',14,'TitleFontSizeMultiplier',1.3)
end

annotation(gcf,'textbox',...
    [0.08 0.89 0.015 0.05],...
    'String','a',...
    'FontSize',24,...
    'FitBoxToText','on',...
    'EdgeColor','none');

annotation(gcf,'textbox',...
    [0.08 0.55 0.02 0.075],...
    'String',{'b'},...
    'FontSize',24,...
    'FitBoxToText','on',...
    'EdgeColor','none');

annotation(gcf,'textbox',...
    [0.08 0.25 0.025 0.0725],...
    'String','c',...
    'FontSize',24,...
    'FitBoxToText','on',...
    'EdgeColor','none');


savefig('Fig1.fig')
close

%% Figure 1b. Average TTV, ERSP, ITC time courses for each band
figure
for c = 1:settings.nfreqs
    subplot(2,3,c)
    t = linspace(0,length(settings.real.poststim)/settings.srate,length(settings.real.poststim));
    plot(t,NormOntoRange(mean(mean(allmeas{c}.ttv.real,3),1),[0 1]),'b','LineWidth',1.5)
    hold on
    plot(t,NormOntoRange(mean(mean(allmeas{c}.ersp.real,3),1),[0 1]),'r','LineWidth',1.5)
    plot(t,1-NormOntoRange(mean(mean(allmeas{c}.itc.real,3),1),[0 1]),'g','LineWidth',1.5)
    
    title(fbands{c})
    legend({'TTV','ERSP','1/ITC'})
    xlabel('Time (s)')
    ylabel('Normalized Units')
end
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)/2 pos(3)*2 pos(4)*1.5]) %make figure larger

ax = findall(gcf,'Type','Axes');
for c = 1:length(ax)
    set(ax(c),'FontSize',14,'TitleFontSizeMultiplier',1.3)
end

savefig('Fig1b.fig')
close

%% Figure 1c. Time course similarity topoplots

figure

for c = 1:6
    subplot(2,3,c)
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,(alloutputs.itc.distrealreal-alloutputs.ttv.diststats{c}.mask > 0),settings.datasetinfo.label);
    else
        topoplot((alloutputs.itc.distrealreal-alloutputs.ttv.diststats{c}.mask > 0),settings.layout);
    end
    title(fbands{c})
    set(gca,'TitleFontSizeMultiplier',1.3)
end
savefig('Fig1c.fig')
close


%% Figure 2a. Correlation between TTVindex and ERSPindex, ITCindex

%designed for 6 frequency bands

figure
subplot(4,7,[8 9 15 16])
t = linspace(0,length(settings.real.poststim)/settings.srate,length(settings.real.poststim));
x2 = [t(settings.aucindex), fliplr(t(settings.aucindex))];
inBetween = [mean(mean(allmeas{1}.ttv.real(:,settings.aucindex,:),3),1), ...
    fliplr(mean(mean(allmeas{1}.ttv.pseudo(:,settings.aucindex,:),3),1))];
fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(t,mean(mean(allmeas{1}.ttv.real,3),1),'b','LineWidth',1.5); %assumes 1 is broadband
plot(t,mean(mean(allmeas{1}.ttv.pseudo,3),1),'b--','LineWidth',1.5)
xlabel('Time (s)')
ylabelunits(settings)
p = get(gca,'position');
set(gca,'FontSize',14,'position',p+[-0.02 0 -0.04 0],'TitleFontSizeMultiplier',1.3)


subplot(4,7,[3 4 10 11])
t = linspace(0,length(settings.real.poststim)/settings.srate,length(settings.real.poststim));
x2 = [t(settings.aucindex), fliplr(t(settings.aucindex))];
inBetween = [mean(mean(allmeas{1}.ersp.real(:,settings.aucindex,:),3),1), ...
    fliplr(mean(mean(allmeas{1}.ersp.pseudo(:,settings.aucindex,:),3),1))];
fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(t,mean(mean(allmeas{1}.ersp.real,3),1),'r','LineWidth',1.5); %assumes 1 is broadband
plot(t,mean(mean(allmeas{1}.ersp.pseudo,3),1),'r--','LineWidth',1.5)

xlabel('Time (s)')
ylabelunits(settings)
p = get(gca,'position');
set(gca,'FontSize',14,'position',p+[0 0 -0.02 0],'TitleFontSizeMultiplier',1.3)

subplot(4,7,[17 18 24 25])
t = linspace(0,length(settings.real.poststim)/settings.srate,length(settings.real.poststim));
x2 = [t(settings.aucindex), fliplr(t(settings.aucindex))];
inBetween = [mean(mean(allmeas{1}.itc.real(:,settings.aucindex,:),3),1), ...
    fliplr(mean(mean(allmeas{1}.itc.pseudo(:,settings.aucindex,:),3),1))];
fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(t,mean(mean(allmeas{1}.itc.real,3),1),'g','LineWidth',1.5); %assumes 1 is broadband
plot(t,mean(mean(allmeas{1}.itc.pseudo,3),1),'g--','LineWidth',1.5)

xlabel('Time (s)')
ylabelunits(settings)
p = get(gca,'position');
set(gca,'FontSize',14,'position',p+[0 0 -0.02 0],'TitleFontSizeMultiplier',1.3)

%% Figure 2c. Topoplots of index correlation
plotindex1 = [5 6 7 12 13 14];
plotindex2 = [19 20 21 26 27 28];
for c = 1:settings.nfreqs
    subplot(4,7,plotindex1(c))
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,alloutputs.ersp.indexr(:,c).*(alloutputs.ersp.corrstats{c}.mask > 0),settings.datasetinfo.label);
    else
        topoplot(alloutputs.ersp.indexr(:,c).*(alloutputs.ersp.corrstats{c}.mask > 0),settings.layout);
    end
    title(fbands{c})
    
    subplot(4,7,plotindex2(c))
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,alloutputs.itc.indexr(:,c).*(alloutputs.itc.corrstats{c}.mask > 0),settings.datasetinfo.label);
    else
        topoplot(alloutputs.itc.indexr(:,c).*(alloutputs.itc.corrstats{c}.mask > 0),settings.layout);
    end
    title(fbands{c})
end
Normalize_Clim(gcf)
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)/2 pos(3)*3 pos(4)*2]) %make figure larger
% Create colorbar
cbar = colorbar('peer',gca,'Position',...
    [0.93056816119078 0.345454545454545 0.0141586494826639 0.353896101065159],'FontSize',12);
cbar.Label.String = 'Spearman correlation';
cbar.Label.FontSize = 14;

% Create arrow
annotation(gcf,'arrow',[0.2775 0.325],...
    [0.5775 0.7225],'LineWidth',2);

% Create arrow
annotation(gcf,'arrow',[0.2775 0.325],...
    [0.4675 0.315],'LineWidth',2);

% Create textbox
annotation(gcf,'textbox',[0.065 0.67 0.015 0.045],'String','a',...
    'FontSize',24,'FitBoxToText','off','EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',[0.31 0.8625 0.0180 0.06285],'String',{'b'},...
    'FontSize',24,'FitBoxToText','off','EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',[0.56 0.8625 0.018 0.06285],'String','c',...
    'FontSize',24,'FitBoxToText','off','EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',[0.56 0.445 0.0180 0.06285],'String','d',...
    'FontSize',24,'FitBoxToText','off','EdgeColor','none');


% Create textbox
annotation(gcf,'textbox',[0.478690 0.889 0.072619 0.03809525],'String','ERSP Index',...
    'FontSize',18,'FitBoxToText','on','EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',[0.2141 0.666666666666667 0.05792 0.0369],...
    'String','TTV Index','FontSize',18,'FitBoxToText','off','EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',[0.48348 0.44894 0.05792 0.03690],'String','ITC Index',...
    'FontSize',18,'FitBoxToText','off','EdgeColor','none');

savefig('Fig2b.fig')
close

%% Figure 2c. Euclidean distance?

%% Figure 4a. NA plots

subplot(4,5,[1 2 6 7])
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,1,:),4),1),'b','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,1,:),4),1),'b--','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,2,:),4),1),'r--','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,2,:),4),1),'r','LineWidth',1.5)
legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
title('TTV')


subplot(4,5,[11 12 16 17])
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,1,:),4),1)...
    -mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,1,:),4),1),'b','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,2,:),4),1)-...
    mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,2,:),4),1),'r','LineWidth',1.5)
legend({'Corrected prestim low','Corrected prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
title('Corrected TTV')

splitmethod = {'amp','cos'};
plotindex{1} = [3 4 5 8 9 10];
plotindex{2} = [13 14 15 18 19 20];
for q = 1:2
    for c = 1:settings.nfreqs
        subplot(4,5,plotindex{q}(c))
        t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
        hold on
        plot(t,mean(mean(allmeas{c}.naddttv.(splitmethod{q}).real(:,:,1,:),4),1)-...
            mean(mean(allmeas{c}.naddttv.(splitmethod{q}).pseudo(:,:,1,:),4),1),'b','LineWidth',1.5)
        
        plot(t,mean(mean(allmeas{c}.naddttv.(splitmethod{q}).real(:,:,2,:),4),1)-...
            mean(mean(allmeas{c}.naddttv.(splitmethod{q}).pseudo(:,:,2,:),4),1),'r','LineWidth',1.5)
        title(fbands{c})
        %legend({'Corrected prestim low','Corrected prestim high'})
        xlabel('Time (s)')
        ylabelunits(settings)
    end
end
savefig('Fig_4a_1.fig')

% ERSP
figure
subplot(4,5,[1 2 6 7])
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
plot(t,mean(mean(allmeas{1}.naddersp.amp.real(:,:,1,:),4),1),'b','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddersp.amp.pseudo(:,:,1,:),4),1),'b--','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddersp.amp.pseudo(:,:,2,:),4),1),'r--','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddersp.amp.real(:,:,2,:),4),1),'r','LineWidth',1.5)
legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
title('Hilbert amplitude')


subplot(4,5,[11 12 16 17])
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
plot(t,mean(mean(allmeas{1}.naddersp.amp.real(:,:,1,:),4),1)...
    -mean(mean(allmeas{1}.naddersp.amp.pseudo(:,:,1,:),4),1),'b','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddersp.amp.real(:,:,2,:),4),1)-...
    mean(mean(allmeas{1}.naddersp.amp.pseudo(:,:,2,:),4),1),'r','LineWidth',1.5)
legend({'Corrected prestim low','Corrected prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
title('Corrected Hilbert amplitude')

splitmethod = {'amp','cos'};
plotindex{1} = [3 4 5 8 9 10];
plotindex{2} = [13 14 15 18 19 20];
for q = 1:2
    for c = 1:settings.nfreqs
        subplot(4,5,plotindex{q}(c))
        t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
        hold on
        plot(t,mean(mean(allmeas{c}.naddersp.(splitmethod{q}).real(:,:,1,:),4),1)-...
            mean(mean(allmeas{c}.naddersp.(splitmethod{q}).pseudo(:,:,1,:),4),1),'b','LineWidth',1.5)
        
        plot(t,mean(mean(allmeas{c}.naddersp.(splitmethod{q}).real(:,:,2,:),4),1)-...
            mean(mean(allmeas{c}.naddersp.(splitmethod{q}).pseudo(:,:,2,:),4),1),'r','LineWidth',1.5)
        title(fbands{c})
        %legend({'Corrected prestim low','Corrected prestim high'})
        xlabel('Time (s)')
        ylabelunits(settings)
    end
end
savefig('Fig_4a_2.fig')

%% Figure 4b. Significance of NAindex for ERSP and TTV

for q = 1:2
    figure
    for c = 1:settings.nfreqs
        ax(c) = subplot(4,3,c);
        if strcmpi(settings.datatype,'MEG')
            ft_topoplot_vec(settings.layout,median(allmeas{c}.nattvindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ttv.stats{c}.mask),settings.datasetinfo.label);
        else
            topoplot(median(allmeas{c}.nattvindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ttv.stats{c}.mask),settings.layout);
        end
        title(fbands{c})
        set(gca,'TitleFontSizeMultiplier',1.3)
        
        ax(c+settings.nfreqs) = subplot(4,3,c+settings.nfreqs);
        if strcmpi(settings.datatype,'MEG')
            ft_topoplot_vec(settings.layout,median(allmeas{c}.naerspindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ersp.stats{c}.mask),settings.datasetinfo.label);
        else
            topoplot(median(allmeas{c}.naerspindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ersp.stats{c}.mask),settings.layout);
        end
        title(fbands{c})
        set(gca,'TitleFontSizeMultiplier',1.3)
    end
    Normalize_Clim(gcf)
    
    %Aesthetics
    cbar = colorbar('peer',ax(1),'Position',...
        [0.865 0.33404 0.01570685 0.38479],'FontSize',12);
    cbar.Label.String = 'Nonadditivity Index';
    cbar.Label.FontSize = 14;
    
    % Create textbox
    annotation(gcf,'textbox',...
        [0.16178 0.861554 0.01947 0.07283],...
        'String','a',...
        'FontSize',24,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    
    % Create textbox
    annotation(gcf,'textbox',...
        [0.16178 0.423188 0.01947 0.07283],...
        'String','b',...
        'FontSize',24,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    
    savefig(['Fig4b_' num2str(q)])
    close
end



%% Figure 5a: Resting state correlations with ERSPindex

figure
if isfield(settings,'rest')
    subplot(2,5,[1 2 6 7])
    nicecorrplot(squeeze(mean(restmeas.rel_bp.vals(4,:,:),2)),squeeze(mean(allmeas{4}.erspindex,1)),{'Relative resting-state band-limited power','Poststimulus ERSPindex'});
    plotindex = [3 4 5 8 9];
    for c = 2:settings.nfreqs
        ax(c) = subplot(2,5,plotindex(c-1));
        if strcmpi(settings.datatype,'MEG')
            ft_topoplot_vec(settings.layout,restmeas.rel_bp.index.r(:,c).*...
                (restmeas.rel_bp.index.stats{c}.mask),settings.datasetinfo.label);
        else
            topoplot(restmeas.rel_bp.index.r(:,c).*...
                (restmeas.rel_bp.index.stats{c}.mask),settings.layout);
        end
        title(settings.tfparams.fbandnames{c})
        set(gca,'TitleFontSizeMultiplier',1.3)
    end
    
    pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)/2 pos(3)*2 pos(4)*1.5])
    % Create colorbar
cbar = colorbar('peer',ax(2),'Position',...
    [0.801837 0.148936 0.020385668 0.268461],'FontSize',12);
cbar.Label.String = 'Spearman Correlation'
cbar.Label.FontSize = 14;

% Create textbox
annotation(gcf,'textbox',...
    [0.0520725 0.870 0.01824 0.05816],...
    'String','a',...
    'FontSize',24,...
    'FitBoxToText','on',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.43410 0.870 0.01824 0.0581
    



],...
    'String','b',...
    'FontSize',24,...
    'FitBoxToText','on',...
    'EdgeColor','none');
    savefig(['Fig_5a.fig'])
    
    %% Figure 5b: Resting state mediation model
    
    %some schematic plot
    % foi = find(~isempty(restmeas.bp.mediation.stats));
    %
    % for q = 1:length(foi)
    %     subplot(1,
    % end
    
    
    
end

FigsToPNG(pwd)

end

function ylabelunits(settings)
switch settings.units
    case 'prcchange'
        ylabel('% change from prestim')
    case 'log'
        ylabel('10*log10 unit change')
    case 'raw'
        ylabel('Change from prestim (respective units)')
    case 'zscore'
        ylabel('Normalized change from prestim')
end
end

