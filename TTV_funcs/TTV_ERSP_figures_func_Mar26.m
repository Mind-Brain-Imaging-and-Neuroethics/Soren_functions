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

%% First make raw figures

%% Figure 1a. Median split - schematic
if strcmpi(settings.datatype,'MEG')
    plotsensor = find(strcmpi('A1',settings.datasetinfo.label));
else
    plotsensor = find(strcmpi('Oz',settings.datasetinfo.label)); %just for Ivar data - where NA is strongest for schematic
end

trange = (prestim_pseudo(1)-settings.srate/5):(poststim_real(end));

%prestimrange = trange(1:(length(trange)-length(poststim_real)-length(prestim_real)));
%poststimrange = trange(201:end);
warning('Hard-coded')
t = linspace(-(poststim_real(1)-trange(1))*(1/settings.srate),length(poststim_real)*(1/settings.srate),length(trange));

subplot(2,1,1)
plot(t,mean(mean(allmeas{4}.naddersp.amp.raw.pseudo(plotsensor,trange,1,:),4),1),...
    'Color',[0.5 0.5 1],'LineWidth',0.5,'HandleVisibility','off')
hold on
plot(t,mean(mean(allmeas{4}.naddersp.amp.raw.pseudo(plotsensor,trange,2,:),4),1),...
    'Color',[1 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off')
plot(t,abs(hilbert(mean(allmeas{4}.naddersp.amp.raw.pseudo(plotsensor,trange,1,:),4))),'b--','LineWidth',2);
plot(t,abs(hilbert(mean(allmeas{4}.naddersp.amp.raw.pseudo(plotsensor,trange,2,:),4))),'r--','LineWidth',2);
FixAxes(gca,14)
ylim = get(gca,'YLim');
line([0 0],ylim,'Color','k','LineWidth',2,'HandleVisibility','off')
line([-1 -1],ylim,'Color','k','LineWidth',2,'LineStyle','--','HandleVisibility','off')
patch([-1.2 -1 -1 -1.2],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
patch([-0.2 0 0 -0.2],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
set(gca,'YLim',ylim)
xlabel('Time (s)')
ylabel('Voltage (uV)')
legend({'Pseudotrial prestim low','Pseudotrial prestim high'},'EdgeColor','none')

warning('Hard-coded')

subplot(2,1,2)
plot(t,mean(mean(allmeas{4}.naddersp.amp.raw.real(plotsensor,trange,1,:),4),1),...
    'Color',[0.5 0.5 1],'LineWidth',0.5,'HandleVisibility','off')
hold on
plot(t,mean(mean(allmeas{4}.naddersp.amp.raw.real(plotsensor,trange,2,:),4),1),...
    'Color',[1 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off')
plot(t,abs(hilbert(mean(allmeas{4}.naddersp.amp.raw.real(plotsensor,trange,1,:),4))),'b','LineWidth',2);
plot(t,abs(hilbert(mean(allmeas{4}.naddersp.amp.raw.real(plotsensor,trange,2,:),4))),'r','LineWidth',2);
FixAxes(gca)
ylim = get(gca,'YLim');
line([0 0],ylim,'Color','k','LineWidth',2,'HandleVisibility','off')
line([-1 -1],ylim,'Color','k','LineWidth',2,'LineStyle','--','HandleVisibility','off')
patch([-1.2 -1 -1 -1.2],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
patch([-0.2 0 0 -0.2],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
set(gca,'YLim',ylim)
%warning('Hard-coded')
xlabel('Time (s)')
ylabel('Voltage (uV)')
legend({'Real trial prestim low','Real trial prestim high'},'EdgeColor','none')

savefig('Fig1a.fig')
close

%% Figure 1b: Normalized pseudotrial and real trial ERSP

figure
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
stdshade(t,squeeze(mean(allmeas{4}.naddersp.amp.real(:,:,1,:),1)),'b',0,1,'sem')
stdshade(t,squeeze(mean(allmeas{4}.naddersp.amp.pseudo(:,:,1,:),1)),'b--',0,1,'sem')
stdshade(t,squeeze(mean(allmeas{4}.naddersp.amp.real(:,:,2,:),1)),'r',0,1,'sem');
stdshade(t,squeeze(mean(allmeas{4}.naddersp.amp.pseudo(:,:,2,:),1)),'r--',0,1,'sem')
%FillBetween(t,mean(mean(allmeas{1}.nadderp.real(:,:,1,:),4),1)-...
%    mean(mean(allmeas{1}.nadderp.pseudo(:,:,1,:),4),1),mean(mean(allmeas{1}.nadderp.real(:,:,2,:),4),1)-...
%    mean(mean(allmeas{1}.nadderp.pseudo(:,:,2,:),4),1));
legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'},'EdgeColor','none')
xlabel('Time (s)')
ylabelunits(settings)
FixAxes(gca,20)

savefig('Fig1b.fig')
close

%% Figure 1c: Corrected real-trial time courses and illustration of pseudotrial index

t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
stdshade(t,squeeze(mean(allmeas{4}.naddersp.amp.real(:,:,1,:),1))-...
    squeeze(mean(allmeas{4}.naddersp.amp.pseudo(:,:,1,:),1)),'b',0,1,'sem')
stdshade(t,squeeze(mean(allmeas{4}.naddersp.amp.real(:,:,2,:),1))-...
    squeeze(mean(allmeas{4}.naddersp.amp.pseudo(:,:,2,:),1)),'r',0,1,'sem')
FillBetween(t,mean(mean(allmeas{4}.naddersp.amp.real(:,:,1,:),4),1)-...
    mean(mean(allmeas{4}.naddersp.amp.pseudo(:,:,1,:),4),1),mean(mean(allmeas{4}.naddersp.amp.real(:,:,2,:),4),1)-...
    mean(mean(allmeas{4}.naddersp.amp.pseudo(:,:,2,:),4),1));
legend({'Corrected prestim low','Corrected prestim high'},'EdgeColor','none')
xlabel('Time (s)')
ylabelunits(settings)
FixAxes(gca,20)

savefig('Fig1c.fig')
close

%% Figure 1d: Nonadditivity of ERSP in different frequency bands
for c = 1:settings.nfreqs
    subplot(1,6,c)
    t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
    hold on
    stdshade(t,squeeze(mean(allmeas{c}.naddersp.amp.real(:,:,1,:),1)-...
        mean(allmeas{c}.naddersp.amp.pseudo(:,:,1,:),1)),'b',0.15,1,'sem');
    
    stdshade(t,squeeze(mean(allmeas{c}.naddersp.amp.real(:,:,2,:),1)-...
        mean(allmeas{c}.naddersp.amp.pseudo(:,:,2,:),1)),'r',0.15,1,'sem');
    title(fbands{c})
    %legend({'Corrected prestim low','Corrected prestim high'})
    xlabel('Time (s)')
    ylabelunits(settings)
    FixAxes(gca)
    set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
end
savefig('Fig1d.fig')
close

%% Figure 1e. Significance of NAindex for ERSP

figure
for c = 1:settings.nfreqs
    ax(c) = subplot(1,6,c);
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,median(allmeas{c}.naerspindex.amp,2).*...
            (alloutputs.([splitmethod{q} '_na']).ersp.stats{c}.mask),settings.datasetinfo.label);
    else
        cluster_topoplot(median(allmeas{c}.naerspindex.amp,2),...
            settings.layout,alloutputs.amp_na.ersp.sig(c,:)',(alloutputs.amp_na.ersp.stats{c}.mask));
    end
    colormap(lkcmap2)
    title(fbands{c})
    set(gca,'TitleFontSizeMultiplier',1.3)
end
Normalize_Clim(gcf)

%Aesthetics
cbar = colorbar('peer',ax(1),'FontSize',12);
cbar.Label.String = 'Pseudotrial index';
cbar.Label.FontSize = 14;

savefig('Fig1e.fig')
close
    
%% Figure 5a: TTV of ERSP

figure
t = linspace(-length(settings.real.prestim)*(1/settings.srate),length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim)+length(settings.real.prestim));
prestimdata = 100*(allmeas{4}.raw.ttversp(:,prestim_real,:)-mean(allmeas{4}.raw.ttversp(:,prestim_real,:),2))./...
    mean(allmeas{4}.raw.ttversp(:,prestim_real,:),2); %assuming percent change units
prestimdata = mean(mean(prestimdata,3),1);
poststimdata = 100*(allmeas{4}.ttversp.real)./...
    mean(allmeas{4}.raw.ttversp(:,prestim_real,:),2); %REMOVE THIS IF YOU RECALCULATE
poststimdata = mean(mean(poststimdata,3),1);
plotdata = [prestimdata poststimdata];
FillBetween(t((length(settings.real.prestim)+1):end),poststimdata,...
    zeros(1,length(poststimdata)));
hold on
plot(t,plotdata,'b','LineWidth',3)
plot(t,zeros(1,length(plotdata)),'k--','LineWidth',1.5)
xlabel('Time (s)')
ylabel('% change of TTV of ERSP')
ylim = get(gca,'YLim');
line([0 0],ylim,'Color',[0.5 0.5 0.5],'LineWidth',2)
set(gca,'YLim',ylim)
FixAxes(gca)
set(gca,'FontSize',20)

aucdata = 100*(allmeas{4}.ttversp.real)./...
    mean(allmeas{4}.raw.ttversp(:,prestim_real,:),2); %REMOVE THIS IF YOU RECALCULATE
topoplot_data = squeeze(trapz(aucdata(:,settings.aucindex,:),2));
p = signrank_mat(topoplot_data,zeros(size(topoplot_data)),2);
plotstats = EasyClusterCorrect_signrank({topoplot_data,zeros(size(topoplot_data))},settings.datasetinfo);
axes('position',[0.15 0.15 0.25 0.25])
cluster_topoplot(mean(topoplot_data,2),settings.layout,p,plotstats.mask)
colorbar('FontSize',12)
title(['p = ' num2str(round(plotstats.posclusters.prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
colormap(lkcmap2)
savefig('Fig5a.fig')

%% Figure 5b: Correlation of TTVERSP with ERSP NA index

nicecorrplot(mean(allmeas{4}.ttverspindex,1),mean(allmeas{4}.naerspindex.amp,1),{'TTV of ERSP Index','Nonadditivity index'})
FixAxes(gca,20)
axes('position',[0.15 0.15 0.25 0.25])
cluster_topoplot(alloutputs.ttversp.r(:,4),settings.layout,alloutputs.ttversp.p(:,4),alloutputs.ttversp.corrstats{4}.mask);
colorbar('FontSize',12)
title(['p = ' num2str(round(plotstats.posclusters.prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
savefig('Fig5b.fig')

%% Figure 3a. Raw TTV, Amplitude, ITC time courses

figure

allrange = min(settings.pseudo.prestim):max(settings.real.poststim);
onset = min(settings.real.poststim);
for c = 1:settings.nfreqs
    ax(c) = subplot(3,settings.nfreqs,c);
    %t = linspace((1/settings.srate)*(min(allrange)-onset),(1/settings.srate)*(max(allrange)-onset),length(allrange));
    t = linspace(0,(1/settings.srate)*length(poststim_pseudo),length(poststim_pseudo));
    stdshade(t,squeeze(mean(allmeas{c}.raw.sd(:,poststim_real,:),1)),'b',0.15,1);
    hold on
    stdshade(t,squeeze(mean(allmeas{c}.raw.sd(:,poststim_pseudo,:),1)),'b--',0.15,1);
    if c == 1
        ylabel('Across-trial SD')
    end
    title(fbands{c})
    FixAxes(gca)
    ax(c).YLabel.FontSize = 18;
    
    ax(c+settings.nfreqs) = subplot(3,settings.nfreqs,c+settings.nfreqs);
    stdshade(t,squeeze(mean(allmeas{c}.raw.ersp(:,poststim_real,:),1)),'r',0.15,1)
    hold on
    stdshade(t,squeeze(mean(allmeas{c}.raw.ersp(:,poststim_pseudo,:),1)),'r--',0.15,1)
    if c == 1
        ylabel('Hilbert amplitude')
    end
    FixAxes(gca)
    ax(c+settings.nfreqs).YLabel.FontSize = 18;

    
    ax(c+2*settings.nfreqs) = subplot(3,settings.nfreqs,c+2*settings.nfreqs);
    stdshade(t,squeeze(mean(allmeas{c}.raw.itc(:,poststim_real,:),1)),'g',0.15,1)
    hold on
    stdshade(t,squeeze(mean(allmeas{c}.raw.itc(:,poststim_pseudo,:),1)),'g--',0.15,1)
    xlabel('Time (s)')
    
    if c == 1
        ylabel('Inter-trial coherence')
    end
    FixAxes(gca)
    ax(c+2*settings.nfreqs).YLabel.FontSize = 18;
    
end
%Aesthetics
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)/2 pos(3)*3 pos(4)*2]) %make figure larger
%AddFigureLabel(ax(1),'a')
%AddFigureLabel(ax(7),'b')
%AddFigureLabel(ax(13),'c')


% ax = findall(gcf,'Type','Axes');
% for c = 1:length(ax)
% FixAxes(ax(c))
% end

savefig('Fig3a.fig')
close

%% Figure 3b. Time course similarity boxplots and topoplots

figure
colors = get(groot,'defaultAxesColorOrder');
newcolors = colors(1:6,:);
set(groot,'defaultAxesColorOrder',newcolors)
violinplot([squeeze(mean(alloutputs.ersp.distrealreal(:,:,:),1)) squeeze(mean(alloutputs.itc.distrealreal(:,:,:),1))],[])
set(groot,'defaultAxesColorOrder',colors);
set(gca,'XTickLabel',{'ERSP','ITC'},'XTick',[3.5 9.5])
FixAxes(gca,20)
ylabel('Euclidean Distance')
a = findobj('Parent',gca,'Type','Patch');
legend(a(fliplr(2:2:12)),fbands,'FontSize',12,'EdgeColor','w')
line([6.5 6.5],get(gca,'YLim'),'Color','k','LineWidth',1,'LineStyle','--','HandleVisibility','off')
savefig('Fig2b_1.fig')


figure
for c = 1:6
    subplot(2,3,c)
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,(mean(alloutputs.itc.distrealreal(:,:,c),2)-mean(alloutputs.ersp.distrealreal(:,:,c),2)).*...
            (alloutputs.ttv.diststats{c}.mask > 0),settings.datasetinfo.label);
    else
        cluster_topoplot((mean(alloutputs.itc.distrealreal(:,:,c),2)-mean(alloutputs.ersp.distrealreal(:,:,c),2)),...
            settings.layout,alloutputs.ttv.sigerspvitc(c,:),alloutputs.ttv.diststats{c}.mask);
    end
    colormap(lkcmap2)
    title(fbands{c})
    set(gca,'TitleFontSizeMultiplier',1.3)
end
cbar = colorbar('peer',gca,'Position',...
    [0.90287 0.37857 0.02429 0.3071429],'FontSize',12);
cbar.Label.String = 'Euclidean Distance Difference';
cbar.Label.FontSize = 14;
savefig('Fig2b_2.fig')
close


%% Figure 3c. Scatter plot of index correlation for ERSP, ITC broadband

figure
subplot(1,2,1)
nicecorrplot(mean(allmeas{1}.ttvindex(:,:),1),mean(allmeas{1}.erspindex(:,:),1),{'TTV Index','ERSP Index'});
FixAxes(gca)
subplot(1,2,2)
nicecorrplot(mean(allmeas{1}.ttvindex(:,:),1),mean(allmeas{1}.itcindex(:,:),1),{'TTV Index','ITC Index'});
FixAxes(gca)
savefig('Fig3c.fig')

%% Figure 4a & 4b: Results of median split ERP - broadband only
%plotindex{1} = [3 4 5 8 9 10];
%plotindex{2} = [13 14 15 18 19 20];
figure
subplot(1,2,1)
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
stdshade(t,squeeze(mean(allmeas{1}.nadderp.real(:,:,1,:),1)),'b',0,1,'sem')
stdshade(t,squeeze(mean(allmeas{1}.nadderp.pseudo(:,:,1,:),1)),'b--',0,1,'sem')
stdshade(t,squeeze(mean(allmeas{1}.nadderp.real(:,:,2,:),1)),'r',0,1,'sem');
stdshade(t,squeeze(mean(allmeas{1}.nadderp.pseudo(:,:,2,:),1)),'r--',0,1,'sem')
%FillBetween(t,mean(mean(allmeas{1}.nadderp.real(:,:,1,:),4),1)-...
%    mean(mean(allmeas{1}.nadderp.pseudo(:,:,1,:),4),1),mean(mean(allmeas{1}.nadderp.real(:,:,2,:),4),1)-...
%    mean(mean(allmeas{1}.nadderp.pseudo(:,:,2,:),4),1));
legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
FixAxes(gca,20)

subplot(1,2,2)
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
stdshade(t,squeeze(mean(allmeas{1}.nadderp.real(:,:,1,:),1))-...
    squeeze(mean(allmeas{1}.nadderp.pseudo(:,:,1,:),1)),'b',0,1,'sem')
stdshade(t,squeeze(mean(allmeas{1}.nadderp.real(:,:,2,:),1))-...
    squeeze(mean(allmeas{1}.nadderp.pseudo(:,:,2,:),1)),'r',0,1,'sem')
FillBetween(t,mean(mean(allmeas{1}.nadderp.real(:,:,1,:),4),1)-...
    mean(mean(allmeas{1}.nadderp.pseudo(:,:,1,:),4),1),mean(mean(allmeas{1}.nadderp.real(:,:,2,:),4),1)-...
    mean(mean(allmeas{1}.nadderp.pseudo(:,:,2,:),4),1));
legend({'Corrected prestim low','Corrected prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
FixAxes(gca,20)
axes('position',[0.75 0.135 0.15 0.2])
cluster_topoplot(mean(allmeas{1}.naerpindex,2),settings.layout,alloutputs.erp_na.p(:,1),alloutputs.erp.stats{1}.mask)
title('No significant clusters','FontWeight','normal','FontSize',14)
colormap(lkcmap)
colorbar('WestOutside')
%set(gcf,'Color','w')
%set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
savefig(['Fig2b.fig'])
close

%% Figure 5a: Poststim TTV decrease

figure
t = linspace(-length(settings.real.prestim)*(1/settings.srate),length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim)+length(settings.real.prestim));
prestimdata = 100*(allmeas{1}.raw.sd(:,prestim_real,:)-mean(allmeas{1}.raw.sd(:,prestim_real,:),2))./...
    mean(allmeas{1}.raw.sd(:,prestim_real,:),2);
prestimdata = mean(mean(prestimdata,3),1);
poststimdata = mean(mean(allmeas{1}.ttv.real,3),1)
plotdata = [prestimdata poststimdata];
FillBetween(t((length(prestimdata)+1):end),poststimdata,zeros(1,length(poststimdata)));
hold on
plot(t,plotdata,'b','LineWidth',1.5)
plot(t,zeros(1,length(plotdata)),'k--','LineWidth',1.5)
xlabel('Time (s)')
ylabel('% change of Trial-to-Trial SD')
 ylim = get(gca,'YLim');
line([0 0],ylim,'Color',[0.5 0.5 0.5],'LineWidth',2)
set(gca,'YLim',ylim)
FixAxes(gca)
set(gca,'FontSize',20)

topoplot_data = squeeze(trapz(allmeas{1}.ttv.real(:,settings.aucindex,:),2));
p = signrank_mat(topoplot_data,zeros(size(topoplot_data)),2);
%plotstats = EasyClusterCorrect_signrank({topoplot_data,zeros(size(topoplot_data))},settings.datasetinfo);
axes('position',[0.15 0.15 0.25 0.25])
cluster_topoplot(mean(topoplot_data,2),settings.layout,p,plotstats.mask,2)
%topoplot(mean(topoplot_data,2),settings.layout,'emarker2',{find(plotstats.mask),'o','w',3,1})
colorbar('FontSize',12)
title(['p = ' num2str(round(plotstats.posclusters.prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
colormap(lkcmap)
savefig('Fig2a.fig')

%% Figure 5b. Correlation of NAindex with TTVindex

figure
nicecorrplot(mean(allmeas{1}.naerpindex,1)',mean(allmeas{1}.ttvindex,1)',{'Nonadditivity index','TTV index'})
axes('position',[0.18 0.15 0.25 0.25])
cluster_topoplot(alloutputs.erp_na.r(:,1),settings.layout,alloutputs.erp_na.p(:,1),alloutputs.erp.corrstats{1}.mask)
colorbar('EastOutside','FontSize',12)
colormap(lkcmap)
title('No significant clusters','FontSize',14,'FontWeight','normal')


%% Figure 6a: Topoplots of resting state correlations with ERSPindex

figure

for c = 2:6
    ax(c) = subplot(2,3,c-1);
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,restmeas.rel_bp.index.r.subject(:,c).*...
            (restmeas.rel_bp.index.stats{c}.mask),settings.datasetinfo.label);
    else
        cluster_topoplot(restmeas.rel_bp.index.r.subject(:,c),settings.layout,...
            restmeas.rel_bp.index.p.subject(:,c),(restmeas.rel_bp.index.stats{c}.mask));
    end
    title(settings.tfparams.fbandnames{c})
end
set(gca,'TitleFontSizeMultiplier',1.3)
cbar = colorbar('peer',gca,'FontSize',12);
cbar.Label.String = "Spearman's rho";
cbar.Label.FontSize = 14;
savefig('Fig6a.fig')
close

%% Figure 6b: Electrode-based correlation with ERSPindex
figure
bandindex = find(strcmpi(fbands,'Alpha'));
nicecorrplot(mean(squeeze(restmeas.rel_bp.vals(bandindex,:,:)),1),mean(allmeas{bandindex}.erspindex,1),{'Resting-state relative alpha power','Alpha ERSPindex'})
FixAxes(gca)
savefig('Fig6b.fig')
close

%% Figure 6c: Topoplots of resting state correlations with ERSP NAindex
for c = 2:6
    ax(c) = subplot(2,3,c-1);
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,restmeas.rel_bp.naindex.r.subject(:,c).*...
            (restmeas.rel_bp.naindex.stats{c}.mask),settings.datasetinfo.label);
    else
        cluster_topoplot(restmeas.rel_bp.naindex.r.subject(:,c),...
            settings.layout,restmeas.rel_bp.naindex.p.subject(:,c),(restmeas.rel_bp.naindex.stats{c}.mask));
    end
    title(settings.tfparams.fbandnames{c})
end
set(gca,'TitleFontSizeMultiplier',1.3)
cbar = colorbar('peer',gca','FontSize',12);
cbar.Label.String = "Spearman's rho";
cbar.Label.FontSize = 14;
savefig('Fig6c.fig')
close

%% Figure 6d: Electrode-based correlation with ERSP NAindex

figure
bandindex = find(strcmpi(fbands,'Alpha'));
nicecorrplot(mean(squeeze(restmeas.rel_bp.vals(bandindex,:,:)),1),mean(allmeas{bandindex}.naerspindex.amp,1),{'Resting-state relative alpha power','Alpha Nonadditivity Index'})
FixAxes(gca)
savefig('Fig6d.fig')
close

%% Figure 6e: Mediation model
figure
mediationAnalysis0(double(mean(allmeas{bandindex}.erspindex(find(restmeas.rel_bp.prestim.stats{bandindex}.mask),:),1))',...
    double(squeeze(mean(restmeas.rel_bp.vals(bandindex,find(restmeas.rel_bp.prestim.stats{bandindex}.mask),:),2))),...
    double(mean(restmeas.prestimamp.rel{bandindex}(find(restmeas.rel_bp.prestim.stats{bandindex}.mask),:),1))');
savefig('Fig6e.fig')
close

%% Combining the sub-figures

%% Figure 1

p = panel('no-manage-font')
p.pack('h',{1/3 1/3 1/3})
p(1).pack('v',{50 50})
Fig1a = open('Fig1a.fig')
figaxes = findobj('Parent',Fig1a,'Type','axes');
p(1,1).select(figaxes(2))
FixAxes(p(1,1).axis,16)
p(1,2).select(figaxes(1))
FixAxes(p(1,2).axis,16)
close(Fig1a)

Fig1bc = open('Fig1bc.fig')
figaxes = findobj('Parent',Fig1bc,'Type','axes');
p(2).select(figaxes(3))
FixAxes(p(2).axis,18)
p(3).pack()
p(3).pack({[0.7 0.02 0.3 0.3]})
p(3,1).select(figaxes(2))
FixAxes(p(3,1).axis,18)
p(3,2).select(figaxes(1))
colormap(lkcmap)
set(p(3,2).axis,'FontSize',14)

p.margin = [22 22 10 5]
p.de.margin = [5 5 5 5];
p(2).marginleft = 25;
p(3).marginleft = 25;

%% Figure 2

p = panel('no-manage-font')
p.pack('h',{50 50})
p(1).pack()
p(1).pack({[0.02 0.02 0.3 0.3]})
Fig2a = open('Fig2a.fig')
figaxes = findobj('Parent',Fig2a,'Type','axes');
p(1,1).select(figaxes(2))
FixAxes(p(1,1).axis,18)
p(1,2).select(figaxes(1))
set(p(1,2).axis,'FontSize',14)
colormap(lkcmap)
close(Fig2a)
AddFigureLabel(p(1,1).axis,'A')

p(2).pack()
p(2).pack({[0.02 0.02 0.3 0.3]})
Fig2b = open('Fig2b.fig')
figaxes = findobj('Parent',Fig2b,'Type','axes');
p(2,1).select(figaxes(2))
FixAxes(p(2,1).axis,18)
p(2,2).select(figaxes(1))
set(p(2,2).axis,'FontSize',14)
colormap(lkcmap)

p.margin = [25 20 5 5];
p.de.margin = [5 5 5 5];
p(2).marginleft = 30;




%% Figure 3

figure
p = panel('no-manage-font')
p.pack('v',{40 25 15 20})
p(1).pack(3,6)
Fig3a = open('Fig3a.fig');
figaxes = findobj('Parent',Fig3a,'Type','axes');
figaxes = rot90(rot90(reshape(figaxes,3,6)));
for c = 1:3
    for cc = 1:6
        p(1,c,cc).select(figaxes(c,cc));
        FixAxes(p(1,c,cc).axis,12);
    end
end
close(Fig3a)

Fig3b = open('Fig3b_1.fig')
figaxes = findobj('Parent',Fig3b,'Type','axes');
p(2).select(figaxes)
FixAxes(p(2).axis,18)
close(Fig3b)

Fig3b = open('Fig3b_2.fig')
figaxes = findobj('Parent',Fig3b,'Type','axes');
figaxes = flipud(figaxes);
cbar = findobj('Parent',Fig3b','Type','ColorBar')
p(3).pack('h',{95 3 2})
p(3,1).pack('h',6)
for c = 1:6
    p(3,1,c).select(figaxes(c));
end
p(3,2).select(cbar)
close(Fig3b)

Fig3c = open('Fig3c.fig')
figaxes = findobj('Parent',Fig3c,'Type','axes');
p(4).pack('h',{50 50})
p(4,1).select(figaxes(2))
FixAxes(p(4,1).axis,16)
p(4,2).select(figaxes(1))
FixAxes(p(4,2).axis,16)
p.margin = [28 20 8 10]
p.de.margin = [5 5 5 5];
p(1).margin = [5 15 5 5];
p(1).de.margin = [9 9 5 5];
p(2).marginbottom = 18;
p(4).de.marginleft = 20;


%% Figure 4

figure 
p = panel('no-manage-font');
p.pack('h',{50 50})
p(1).pack('v',{1/3 1/3 1/3});
p(2).pack('v',{1/3 1/3 1/3});
Fig4a = open('Fig4a.fig')
figaxes = findobj('Parent',Fig4a,'Type','Axes');
p(1,1).select(figaxes(2));
AddFigureLabel(p(1,1).axis,'A')
p(2,1).select(figaxes(1));
AddFigureLabel(p(1,1).axis,'B')
close(Fig4a)

for c = 1:2
   Fig4b = open(['Fig4b_' num2str(c) '.fig']);
   figaxes = findobj('Parent',Fig4b,'Type','Axes');
   figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
   p(c,2).pack(2,3)
   for cc = 1:2
       for ccc = 1:3
            p(c,2,cc,ccc).select(figaxes(cc,ccc));
            FixAxes(p(c,2,cc,ccc).axis,12)
            set(gca,'TitleFontSizeMultiplier',1.1)
            if cc == 1
               set(gca,'XLabel',[]) 
            end
            if ccc ~= 1
               set(gca,'YLabel',[]) 
            end
       end
   end
   
    close(Fig4b)   
end


for c = 1:2
   Fig4c = open(['Fig4c_' num2str(c) '.fig']);
   figaxes = findobj('Parent',Fig4c,'Type','Axes');
   figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
   p(c,3).pack(2,3)
   for cc = 1:2
       for ccc = 1:3
            p(c,3,cc,ccc).select(figaxes(cc,ccc));
            FixAxes(p(c,3,cc,ccc).axis,12)
       end
   end
   close(Fig4c)
end
p.margin = [30 10 5 15];
p.de.margin = [5 5 5 5];
p(1).marginright = 30;
p(1,1).marginbottom = 25;
p(2,1).marginbottom = 25;
p(1,2).marginbottom = 25;
p(2,2).marginbottom = 25;
p(1,2).de.margin = [10 15 5 5];
p(2,2).de.margin = [10 15 5 5];
p(1,3).de.margin = [5 10 5 5];
p(2,3).de.margin = [5 10 5 5];
Normalize_Clim(gcf)
colormap(lkcmap)

savefig('Fig4.fig')
p.export('Fig4','-rp')
save('Panel4.mat','p')
clear p
close


%% Figure 5

figure
p = panel('no-manage-font');
p.pack('h',{50 50})
p(1).pack()
p(1).pack({[0.05 0.05 0.3 0.3]})
Fig5a = open('Fig5a.fig')
figaxes = findobj('Parent',Fig5a,'Type','axes')
p(1,1).select(figaxes(2));
FixAxes(p(1,1).axis,18)
p(1,2).select(figaxes(1));
set(p(1,2).axis,'FontSize',14)
close(Fig5a)
AddFigureLabel(p(1,1).axis,'A')

%FixAxes(p(1).axis,18)
Fig5b = open('Fig5b.fig')
figaxes = findobj('Parent',Fig5b,'Type','axes')
p(2).pack()
p(2).pack({[0.05 0.05 0.3 0.3]})
p(2,1).select(figaxes(2));
FixAxes(p(2,1).axis,18)
AddFigureLabel(p(2,1).axis,'B')
p(2,2).select(figaxes(1));
set(p(2,2).axis,'FontSize',14)

p.marginbottom = 25;
p.marginleft = 30;
p(2).marginleft = 35;
p(1,2).margin = [5 5 5 5];

savefig('Fig5.fig')
p.export('Fig5','-rp')
save('Panel5.mat','p')
clear p
close

%% Figure 6

figure
p = panel('no-manage-font');
p.pack('v',{26 26 48})
p(1).pack('h',{1/3 2/3})
p(1,2).pack(2,3)
Fig6a = open('Fig6a.fig');
figaxes = findobj('Parent',Fig6a,'Type','axes');
cbar = findobj('Parent',Fig6a,'Type','ColorBar');
%figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
plotindex = [5 4 3; 2 1 0];
for c = 1:2
    for cc = 1:3
        if c == 2 && cc == 3
            p(1,2,c,cc).pack('h',{30 20 50});
            p(1,2,c,cc,2).select(cbar);
        else
            p(1,2,c,cc).select(figaxes(plotindex(c,cc)));
            FixAxes(p(1,2,c,cc).axis)
        end
    end
end
close(Fig6a)

Fig6b = open('Fig6b.fig');
figaxes = findobj('Parent',Fig6b,'Type','axes');
p(1,1).select(figaxes)
FixAxes(p(1,1).axis)
close(Fig6b)

p(2).pack('h',{1/3 2/3})
p(2,2).pack(2,3)
Fig6c = open('Fig6c.fig');
figaxes = findobj('Parent',Fig6c,'Type','axes');
cbar = findobj('Parent',Fig6c,'Type','ColorBar');
plotindex = [5 4 3; 2 1 0];
%figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
for c = 1:2
    for cc = 1:3
        if c == 2 && cc == 3
            p(2,2,c,cc).pack('h',{30 20 50});
            p(2,2,c,cc,2).select(cbar);
        else
            p(2,2,c,cc).select(figaxes(plotindex(c,cc)));
            FixAxes(p(2,2,c,cc).axis)
        end
    end
end
close(Fig6c)

Fig6d = open('Fig6d.fig');
figaxes = findobj('Parent',Fig6d,'Type','axes');
p(2,1).select(figaxes)
FixAxes(p(2,1).axis)
close(Fig6d)

Fig6e = open('Fig6e.fig');
figaxes = findobj('Parent',Fig6e,'Type','axes');
p(3).select(figaxes)
close(Fig6e)
colormap(lkcmap2)

p.margintop = 10;
p.marginleft = 30;
p.de.margin = [5 5 5 5];
p(1,2).de.margin = [5 10 5 5];
p(1,2).de.marginleft = 5;
p(2).margintop = 20;
p(2).de.margintop = 5;
p(2,2).de.margin = [5 10 5 5];
p(2,2).de.marginleft = 5;
p(3).margintop = 15;

savefig('Fig6.fig')
p.export('Fig6','-rp')
save('Panel6.mat','p')
clear p
close
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

