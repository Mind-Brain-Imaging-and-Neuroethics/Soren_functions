function TTV_ERSP_figures_func(settings)

fbands = settings.tfparams.fbandnames;
load('lkcmap2.mat')

load([settings.outputdir '/' settings.datasetname '_allmeas.mat'])
load([settings.outputdir '/' settings.datasetname '_results.mat'])
if isfield(settings,'rest')
    load([settings.outputdir '/' settings.datasetname '_restmeas.mat'])
end

if strcmpi(settings.datatype,'EEG')
    %eeglab
end

cd([settings.outputdir '/' settings.datasetname '_figures'])

if strcmpi(settings.tfparams.method,'hilbert') || ~isempty(find(contains(settings.steps,'calc')))
    prestim_pseudo = settings.pseudo.prestim;
    prestim_real = settings.real.prestim;
    poststim_pseudo = settings.pseudo.poststim;
    poststim_real = settings.real.poststim;
else
    prestim_pseudo = settings.pseudo.prestim - settings.pseudo.prestim(1)+1+settings.srate/5;
    prestim_real = settings.real.prestim - settings.pseudo.prestim(1)+1+settings.srate/5;
    poststim_pseudo = settings.pseudo.poststim - settings.pseudo.prestim(1)+1+settings.srate/5;
    poststim_real = settings.real.poststim - settings.pseudo.prestim(1)+1+settings.srate/5;
end


%% First make raw figures

%% Figure 1a. Median split - schematic
if strcmpi(settings.datatype,'MEG')
    plotsensor = find(strcmpi('MEG2011',settings.datasetinfo.label));
elseif strcmpi(settings.datatype,'EEG')
    plotsensor = find(strcmpi('Oz',settings.datasetinfo.label)); %just for Ivar data - where NA is strongest for schematic
elseif strcmpi(settings.datatype,'ECoG') && strcmpi(settings.ecog.method,'roi')
    plotsensor = 7;
elseif strcmpi(settings.datatype,'ECoG')
    plotsensor = 1;
end

trange = (prestim_pseudo(1)-settings.srate/5):(poststim_real(end));

%prestimrange = trange(1:(length(trange)-length(poststim_real)-length(prestim_real)));
%poststimrange = trange(201:end);
warning('Hard-coded')
t = linspace(-(poststim_real(1)-trange(1))*(1/settings.srate),length(poststim_real)*(1/settings.srate),length(trange));

plotband = find(strcmpi(fbands,'Alpha'));

subplot(2,1,1)
plot(t,nanmean(nanmean(allmeas{plotband}.naddersp.raw.pseudo(plotsensor,trange,1,:),4),1),...
    'Color',[0.5 0.5 1],'LineWidth',0.5,'HandleVisibility','off')
hold on
plot(t,nanmean(nanmean(allmeas{plotband}.naddersp.raw.pseudo(plotsensor,trange,2,:),4),1),...
    'Color',[1 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off')
plot(t,abs(hilbert(nanmean(allmeas{plotband}.naddersp.raw.pseudo(plotsensor,trange,1,:),4))),'b--','LineWidth',2);
plot(t,abs(hilbert(nanmean(allmeas{plotband}.naddersp.raw.pseudo(plotsensor,trange,2,:),4))),'r--','LineWidth',2);
FixAxes(gca,14)
ylim = get(gca,'YLim');
line([0 0],ylim,'Color','k','LineWidth',2,'HandleVisibility','off')
linepos = t(find(trange == poststim_pseudo(1)));
prestimlen = length(prestim_pseudo)/settings.srate;
line([linepos linepos],ylim,'Color','k','LineWidth',2,'LineStyle','--','HandleVisibility','off')
patch([linepos-prestimlen linepos linepos linepos-prestimlen],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
patch([-prestimlen 0 0 -prestimlen],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
set(gca,'YLim',ylim)
xlabel('Time (s)')
ylabel('Voltage (uV)')
legend({'Pseudotrial prestim low','Pseudotrial prestim high'},'EdgeColor','none')

%warning('Hard-coded')

subplot(2,1,2)
plot(t,nanmean(nanmean(allmeas{plotband}.naddersp.raw.real(plotsensor,trange,1,:),4),1),...
    'Color',[0.5 0.5 1],'LineWidth',0.5,'HandleVisibility','off')
hold on
plot(t,nanmean(nanmean(allmeas{plotband}.naddersp.raw.real(plotsensor,trange,2,:),4),1),...
    'Color',[1 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off')
plot(t,abs(hilbert(nanmean(allmeas{plotband}.naddersp.raw.real(plotsensor,trange,1,:),4))),'b','LineWidth',2);
plot(t,abs(hilbert(nanmean(allmeas{plotband}.naddersp.raw.real(plotsensor,trange,2,:),4))),'r','LineWidth',2);
FixAxes(gca)
ylim = get(gca,'YLim');
line([0 0],ylim,'Color','k','LineWidth',2,'HandleVisibility','off')
line([linepos linepos],ylim,'Color','k','LineWidth',2,'LineStyle','--','HandleVisibility','off')
patch([linepos-prestimlen linepos linepos linepos-prestimlen],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
patch([-prestimlen 0 0 -prestimlen],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
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
stdshade(t,squeeze(nanmean(allmeas{plotband}.naddersp.real(:,:,1,:),1)),'b',0,1,'sem')
stdshade(t,squeeze(nanmean(allmeas{plotband}.naddersp.pseudo(:,:,1,:),1)),'b--',0,1,'sem')
stdshade(t,squeeze(nanmean(allmeas{plotband}.naddersp.real(:,:,2,:),1)),'r',0,1,'sem');
stdshade(t,squeeze(nanmean(allmeas{plotband}.naddersp.pseudo(:,:,2,:),1)),'r--',0,1,'sem')
%FillBetween(t,nanmean(nanmean(allmeas{1}.nadderp.real(:,:,1,:),4),1)-...
%    nanmean(nanmean(allmeas{1}.nadderp.pseudo(:,:,1,:),4),1),nanmean(nanmean(allmeas{1}.nadderp.real(:,:,2,:),4),1)-...
%    nanmean(nanmean(allmeas{1}.nadderp.pseudo(:,:,2,:),4),1));
legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'},'EdgeColor','none')
xlabel('Time (s)')
ylabelunits(settings)
FixAxes(gca,20)

savefig('Fig1b.fig')
close

%% Figure 1c: Corrected real-trial time courses and illustration of pseudotrial index

t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
stdshade(t,squeeze(nanmean(allmeas{plotband}.naddersp.real(:,:,1,:),1))-...
    squeeze(nanmean(allmeas{plotband}.naddersp.pseudo(:,:,1,:),1)),'b',0,1,'sem')
stdshade(t,squeeze(nanmean(allmeas{plotband}.naddersp.real(:,:,2,:),1))-...
    squeeze(nanmean(allmeas{plotband}.naddersp.pseudo(:,:,2,:),1)),'r',0,1,'sem')
FillBetween(t,nanmean(nanmean(allmeas{plotband}.naddersp.real(:,:,1,:),4),1)-...
    nanmean(nanmean(allmeas{plotband}.naddersp.pseudo(:,:,1,:),4),1),nanmean(nanmean(allmeas{plotband}.naddersp.real(:,:,2,:),4),1)-...
    nanmean(nanmean(allmeas{plotband}.naddersp.pseudo(:,:,2,:),4),1));
legend({'Corrected prestim low','Corrected prestim high'},'EdgeColor','none')
xlabel('Time (s)')
ylabelunits(settings)
FixAxes(gca,20)

savefig('Fig1c.fig')
close

%% Figure 1d: Nonadditivity of ERSP in different frequency bands
for c = 1:settings.nfreqs
    subplot(1,settings.nfreqs,c)
    t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
    hold on
    stdshade(t,squeeze(nanmean(allmeas{c}.naddersp.real(:,:,1,:),1)-...
        nanmean(allmeas{c}.naddersp.pseudo(:,:,1,:),1)),'b',0.15,1,'sem');
    
    stdshade(t,squeeze(nanmean(allmeas{c}.naddersp.real(:,:,2,:),1)-...
        nanmean(allmeas{c}.naddersp.pseudo(:,:,2,:),1)),'r',0.15,1,'sem');
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
    ax(c) = subplot(1,settings.nfreqs,c);
    if strcmpi(settings.datatype,'MEG')
        ft_cluster_topoplot(settings.layout,median(allmeas{c}.naerspindex,2),settings.datasetinfo.label,...
            alloutputs.ersp.pt.sig(c,:)',alloutputs.ersp.pt.stats{c}.mask);
    else
        cluster_topoplot(median(allmeas{c}.naerspindex,2),...
            settings.layout,alloutputs.ersp.pt.sig(c,:)',(alloutputs.ersp.pt.stats{c}.mask));
    end
    colormap(lkcmap2)
    title(fbands{c})
    set(gca,'TitleFontSizeMultiplier',1.3)
end
Normalize_Clim(gcf,1)

%Aesthetics
cbar = colorbar('peer',ax(1),'FontSize',12);
cbar.Label.String = 'Pseudotrial index';
cbar.Label.FontSize = 14;

savefig('Fig1e.fig')
close

%% Figure 2a: TTV of ERSP

figure
t = linspace(-length(settings.real.prestim)*(1/settings.srate),length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim)+length(settings.real.prestim));
prestimdata = 100*(allmeas{plotband}.raw.ttversp(:,prestim_real,:)-nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2))./...
    nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2); %assuming percent change units
prestimdata = nanmean(nanmean(prestimdata,3),1);
poststimdata = (allmeas{plotband}.ttversp.real);%./...
%    nanmean(allmeas{4}.raw.ttversp(:,prestim_real,:),2); %REMOVE THIS IF YOU RECALCULATE
poststimdata = nanmean(nanmean(poststimdata,3),1);
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

%aucdata = 100*(allmeas{4}.ttversp.real)./...
%    nanmean(allmeas{4}.raw.ttversp(:,prestim_real,:),2); %REMOVE THIS IF YOU RECALCULATE
topoplot_data = allmeas{plotband}.ttverspindex;
p = alloutputs.ersp.ttv.sig(plotband,:);
plotstats = alloutputs.ersp.ttv.stats{plotband};
%plotstats = EasyClusterCorrect_signrank({topoplot_data,zeros(size(topoplot_data))},settings.datasetinfo);
axes('position',[0.15 0.15 0.25 0.25])
if strcmpi(settings.datatype,'EEG')
cluster_topoplot(nanmean(topoplot_data,2),settings.layout,p,plotstats.mask)
else
   ft_cluster_topoplot(settings.layout,nanmean(topoplot_data,2),settings.datasetinfo.label,p,plotstats.mask); 
end
colorbar('FontSize',12)
if ~isempty(find(plotstats.mask))
    if isfield(plotstats,'posclusters') && ~isempty(plotstats.posclusters)
        title(['p = ' num2str(round(plotstats.posclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
    else
        title(['p = ' num2str(round(plotstats.negclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
    end
end
colormap(lkcmap2)
savefig('Fig2a.fig')
close

%% Figure 2b: Correlation of TTVERSP with ERSP NA index

nicecorrplot(nanmean(allmeas{plotband}.naerspindex,1),nanmean(allmeas{plotband}.ttverspindex,1),{'Pseudotrial-based nonadditivity','TTV-based nonadditivity'})
FixAxes(gca,20)
axes('position',[0.15 0.15 0.25 0.25])
if strcmpi(settings.datatype,'EEG')
cluster_topoplot(alloutputs.ersp.corr.r(:,plotband),settings.layout,alloutputs.ersp.corr.p(:,plotband),alloutputs.ersp.corr.stats{plotband}.mask);
else
   ft_cluster_topoplot(settings.layout,alloutputs.ersp.corr.r(:,plotband),settings.datasetinfo.label,...
       alloutputs.ersp.corr.p(:,plotband),alloutputs.ersp.corr.stats{plotband}.mask);
end
colorbar('FontSize',12)
if ~isempty(find(alloutputs.ersp.corr.stats{plotband}.mask))
    if isfield(alloutputs.ersp.corr.stats{plotband},'posclusters') && ~isempty(alloutputs.ersp.corr.stats{plotband}.posclusters)
        title(['p = ' num2str(round(alloutputs.ersp.corr.stats{plotband}.posclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
    else
        title(['p = ' num2str(round(alloutputs.ersp.corr.stats{plotband}.negclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
        
    end
end
colormap(lkcmap2)
savefig('Fig2b.fig')
close

%% Figure 3a & 3b: Results of median split ERP - broadband only
%plotindex{1} = [3 4 5 8 9 10];
%plotindex{2} = [13 14 15 18 19 20];
figure
p = panel('no-manage-font')
p.pack('h',{50 50});
p(1).select();
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
stdshade(t,squeeze(nanmean(allmeas{1}.nadderp.real(:,:,1,:),1)),'b',0,1,'sem')
stdshade(t,squeeze(nanmean(allmeas{1}.nadderp.pseudo(:,:,1,:),1)),'b--',0,1,'sem')
stdshade(t,squeeze(nanmean(allmeas{1}.nadderp.real(:,:,2,:),1)),'r',0,1,'sem');
stdshade(t,squeeze(nanmean(allmeas{1}.nadderp.pseudo(:,:,2,:),1)),'r--',0,1,'sem')
%FillBetween(t,nanmean(nanmean(allmeas{1}.nadderp.real(:,:,1,:),4),1)-...
%    nanmean(nanmean(allmeas{1}.nadderp.pseudo(:,:,1,:),4),1),nanmean(nanmean(allmeas{1}.nadderp.real(:,:,2,:),4),1)-...
%    nanmean(nanmean(allmeas{1}.nadderp.pseudo(:,:,2,:),4),1));
legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'},'EdgeColor','none')
xlabel('Time (s)')
ylabel('Voltage (uV)')
%ylabelunits(settings)
FixAxes(gca,20)

p(2).pack();
p(2).pack({[0.75 0.02 0.25 0.3]});
p(2,1).select();
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
stdshade(t,squeeze(nanmean(allmeas{1}.nadderp.real(:,:,1,:),1))-...
    squeeze(nanmean(allmeas{1}.nadderp.pseudo(:,:,1,:),1)),'b',0,1,'sem')
stdshade(t,squeeze(nanmean(allmeas{1}.nadderp.real(:,:,2,:),1))-...
    squeeze(nanmean(allmeas{1}.nadderp.pseudo(:,:,2,:),1)),'r',0,1,'sem')
FillBetween(t,nanmean(nanmean(allmeas{1}.nadderp.real(:,:,1,:),4),1)-...
    nanmean(nanmean(allmeas{1}.nadderp.pseudo(:,:,1,:),4),1),nanmean(nanmean(allmeas{1}.nadderp.real(:,:,2,:),4),1)-...
    nanmean(nanmean(allmeas{1}.nadderp.pseudo(:,:,2,:),4),1));
legend({'Corrected prestim low','Corrected prestim high'},'EdgeColor','none')
xlabel('Time (s)')
ylabel('Voltage (uV)')
%ylabelunits(settings)
FixAxes(gca,20)
%axes('position',[0.75 0.135 0.15 0.2])
p(2,2).select();
if strcmpi(settings.datatype,'EEG')
cluster_topoplot(nanmean(allmeas{1}.naerpindex,2),settings.layout,alloutputs.erp.pt.sig(1,:)',alloutputs.erp.pt.stats{1}.mask)
elseif strcmpi(settings.datatype,'MEG')
    ft_cluster_topoplot(settings.layout,nanmean(allmeas{1}.naerpindex,2),settings.datasetinfo.label,...
        alloutputs.erp.pt.sig(1,:)',alloutputs.erp.pt.stats{1}.mask);
end
if ~isempty(find(alloutputs.erp.pt.stats{1}.mask))
    if isfield(alloutputs.erp.pt.stats{1},'posclusters') && ~isempty(alloutputs.erp.pt.stats{1}.posclusters)
        title(['p = ' num2str(round(alloutputs.erp.pt.stats{1}.posclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
    else
        title(['p = ' num2str(round(alloutputs.erp.pt.stats{1}.negclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
    end
end
%title('No significant clusters','FontWeight','normal','FontSize',14)
colormap(lkcmap2)
colorbar('WestOutside')

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*2.5 pos(4)*1.5],'Color','w');

p.margin = [28 22 6 5];
p(1).marginright = 18;
AddFigureLabel(p(1).axis,'A')
AddFigureLabel(p(2,1).axis,'B');


%set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
savefig('Fig3.fig')
export_fig('Fig3.png','-m4')
save('Panel3.mat','p')
close

%% Figure 4a: Poststim TTV decrease

figure
t = linspace(-length(settings.real.prestim)*(1/settings.srate),length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim)+length(settings.real.prestim));
prestimdata = 100*(allmeas{1}.raw.sd(:,prestim_real,:)-nanmean(allmeas{1}.raw.sd(:,prestim_real,:),2))./...
    nanmean(allmeas{1}.raw.sd(:,prestim_real,:),2);
prestimdata = nanmean(nanmean(prestimdata,3),1);
poststimdata = nanmean(nanmean(allmeas{1}.ttv.real,3),1)
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

plotstats = alloutputs.erp.ttv.stats{1};
topoplot_data = squeeze(trapz(allmeas{1}.ttv.real(:,settings.aucindex,:),2));
p = alloutputs.erp.ttv.sig(1,:);
%p = signrank_mat(topoplot_data,zeros(size(topoplot_data)),2);
%plotstats = EasyClusterCorrect_signrank({topoplot_data,zeros(size(topoplot_data))},settings.datasetinfo);
axes('position',[0.15 0.15 0.25 0.25])
if strcmpi(settings.datatype,'EEG')
cluster_topoplot(nanmean(topoplot_data,2),settings.layout,p,plotstats.mask,2);
elseif strcmpi(settings.datatype,'MEG')
   ft_cluster_topoplot(settings.layout,nanmean(topoplot_data,2),settings.datasetinfo.label,p,plotstats.mask)
end
%topoplot(nanmean(topoplot_data,2),settings.layout,'emarker2',{find(plotstats.mask),'o','w',3,1})
colorbar('FontSize',12)
if ~isempty(find(plotstats.mask))
    if isfield(plotstats,'posclusters') && ~isempty(plotstats.posclusters)
        title(['p = ' num2str(round(plotstats.posclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
    else
        title(['p = ' num2str(round(plotstats.negclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
        
    end
end
colormap(lkcmap2)
savefig('Fig4a.fig')
close

%% Figure 4b. Correlation of NAindex with TTVindex

figure
nicecorrplot(nanmean(allmeas{1}.naerpindex,1)',nanmean(allmeas{1}.ttvindex,1)',{'Pseudotrial-based nonadditivity','TTV-based nonadditivity'})
axes('position',[0.18 0.15 0.25 0.25])
if strcmpi(settings.datatype,'EEG')
cluster_topoplot(alloutputs.erp.corr.r(:,1),settings.layout,alloutputs.erp.corr.p(:,1),alloutputs.erp.corr.stats{1}.mask)
elseif strcmpi(settings.datatype,'MEG')
   ft_cluster_topoplot(settings.layout,alloutputs.erp.corr.r(:,1),settings.datasetinfo.label,...
       alloutputs.erp.corr.p(:,1),alloutputs.erp.corr.stats{1}.mask);
end
colorbar('EastOutside','FontSize',12)
colormap(lkcmap2)
plotstats = alloutputs.erp.corr.stats{1};
if ~isempty(find(plotstats.mask))
    if isfield(plotstats,'posclusters') && ~isempty(plotstats.posclusters) && ~isempty(find(extractfield(plotstats.posclusters,'prob') < 0.05))
        title(['p = ' num2str(round(alloutputs.erp.corr.stats{1}.posclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
    else
        title(['p = ' num2str(round(alloutputs.erp.corr.stats{1}.negclusters(1).prob,2,'significant')) '*'],'FontWeight','normal','FontSize',18)
    end
end
savefig('Fig4b.fig')
close


%% Figure 5a. Raw TTV, Amplitude, ITC time courses

figure

allrange = min(settings.pseudo.prestim):max(settings.real.poststim);
onset = min(settings.real.poststim);
for c = 1:settings.nfreqs
    ax(c) = subplot(3,settings.nfreqs,c);
    %t = linspace((1/settings.srate)*(min(allrange)-onset),(1/settings.srate)*(max(allrange)-onset),length(allrange));
    t = linspace(0,(1/settings.srate)*length(poststim_pseudo),length(poststim_pseudo));
    stdshade(t,squeeze(nanmean(allmeas{c}.raw.sd(:,poststim_real,:),1)),'b',0.15,1);
    hold on
    stdshade(t,squeeze(nanmean(allmeas{c}.raw.sd(:,poststim_pseudo,:),1)),'b--',0.15,1);
    if c == 1
        ylabel('Across-trial SD')
    end
    title(fbands{c})
    FixAxes(gca)
    ax(c).YLabel.FontSize = 18;
    
    ax(c+settings.nfreqs) = subplot(3,settings.nfreqs,c+settings.nfreqs);
    stdshade(t,squeeze(nanmean(allmeas{c}.raw.ersp(:,poststim_real,:),1)),'r',0.15,1)
    hold on
    stdshade(t,squeeze(nanmean(allmeas{c}.raw.ersp(:,poststim_pseudo,:),1)),'r--',0.15,1)
    if c == 1
        ylabel('Hilbert amplitude')
    end
    FixAxes(gca)
    ax(c+settings.nfreqs).YLabel.FontSize = 18;
    
    
    ax(c+2*settings.nfreqs) = subplot(3,settings.nfreqs,c+2*settings.nfreqs);
    stdshade(t,squeeze(nanmean(allmeas{c}.raw.itc(:,poststim_real,:),1)),'g',0.15,1)
    hold on
    stdshade(t,squeeze(nanmean(allmeas{c}.raw.itc(:,poststim_pseudo,:),1)),'g--',0.15,1)
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

savefig('Fig5a.fig')
close

%% Figure 5b. Time course similarity boxplots and topoplots

figure
%colors = get(groot,'defaultAxesColorOrder');
colors = jet;
newcolors = colors(round(linspace(1,64,settings.nfreqs)),:);
set(groot,'defaultAxesColorOrder',newcolors)
violinplot([squeeze(nanmean(alloutputs.ersp.distrealreal(:,:,:),1)) squeeze(nanmean(alloutputs.itc.distrealreal(:,:,:),1))],[])
set(groot,'defaultAxesColorOrder',colors);
set(gca,'XTickLabel',{'ERSP','ITC'},'XTick',[settings.nfreqs/2+0.5 3*settings.nfreqs/2+0.5])
FixAxes(gca,20)
ylabel('Euclidean Distance')
a = findobj('Parent',gca,'Type','Patch');
legend(a(fliplr(2:2:(settings.nfreqs*2))),fbands,'FontSize',12,'EdgeColor','w')
line([settings.nfreqs+0.5 settings.nfreqs+0.5],get(gca,'YLim'),'Color','k','LineWidth',1,'LineStyle','--','HandleVisibility','off')
savefig('Fig5b_1.fig')
close

figure
for c = 1:settings.nfreqs
    subplot(2,4,c)
    if strcmpi(settings.datatype,'MEG')
        ft_cluster_topoplot(settings.layout,(nanmean(alloutputs.itc.distrealreal(:,:,c),2)-nanmean(alloutputs.ersp.distrealreal(:,:,c),2)),...
            settings.datasetinfo.label,alloutputs.dist.sigerspvitc(c,:),alloutputs.dist.stats{c}.mask);
    else
        cluster_topoplot((nanmean(alloutputs.itc.distrealreal(:,:,c),2)-nanmean(alloutputs.ersp.distrealreal(:,:,c),2)),...
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
savefig('Fig5b_2.fig')
close


% %% Figure 5c. Scatter plot of index correlation for ERSP, ITC broadband
%
% figure
% subplot(1,2,1)
% nicecorrplot(nanmean(allmeas{1}.ttvindex(:,:),1),nanmean(allmeas{1}.erspindex(:,:),1),{'TTV Index','ERSP Index'});
% FixAxes(gca)
% subplot(1,2,2)
% nicecorrplot(nanmean(allmeas{1}.ttvindex(:,:),1),nanmean(allmeas{1}.itcindex(:,:),1),{'TTV Index','ITC Index'});
% FixAxes(gca)
% savefig('Fig5c.fig')


%% Figure 6a: Topoplots of resting state correlations with ERSPindex
if isfield(settings,'rest')
figure

for c = 2:settings.nfreqs
    ax(c) = subplot(2,settings.nfreqs/2,c-1);
    if strcmpi(settings.datatype,'MEG')
        ft_cluster_topoplot(settings.layout,restmeas.rel_bp.index.r.subject(:,c),settings.datasetinfo.label,...
            restmeas.rel_bp.index.p.subject(:,c),restmeas.rel_bp.index.stats{c}.mask);
    else
        cluster_topoplot(restmeas.rel_bp.index.r.subject(:,c),settings.layout,...
            restmeas.rel_bp.index.p.subject(:,c),(restmeas.rel_bp.index.stats{c}.mask));
    end
    title(settings.tfparams.fbandnames{c})
end
set(gca,'TitleFontSizeMultiplier',1.3)
cbar = colorbar('peer',gca,'FontSize',12);
cbar.Label.String = 'Spearman''s rho';
cbar.Label.FontSize = 14;
colormap(lkcmap2)
newlim = Normalize_Clim(gcf,1);
savefig('Fig6a.fig')
close

%% Figure 6b: Electrode-based correlation with ERSPindex
figure
bandindex = find(strcmpi(fbands,'Alpha'));
nicecorrplot(nanmean(squeeze(restmeas.rel_bp.vals(bandindex,:,:))',1),nanmean(allmeas{bandindex}.erspindex',1),{'Resting-state relative alpha power','Alpha ERSP'})
FixAxes(gca)
savefig('Fig6b.fig')
close

%% Figure 6c: Topoplots of resting state correlations with ERSP NAindex
figure
for c = 2:settings.nfreqs
    ax(c) = subplot(2,settings.nfreqs/2,c-1);
    if strcmpi(settings.datatype,'MEG')
        ft_cluster_topoplot(settings.layout,restmeas.rel_bp.naindex.r.subject(:,c),...
            settings.datasetinfo.label,restmeas.rel_bp.naindex.p.subject(:,c),restmeas.rel_bp.naindex.stats{c}.mask);
    else
        cluster_topoplot(restmeas.rel_bp.naindex.r.subject(:,c),...
            settings.layout,restmeas.rel_bp.naindex.p.subject(:,c),(restmeas.rel_bp.naindex.stats{c}.mask));
    end
    title(settings.tfparams.fbandnames{c})
end
set(gca,'TitleFontSizeMultiplier',1.3)
cbar = colorbar('peer',gca','FontSize',12);
cbar.Label.String = 'Spearman''s rho';
cbar.Label.FontSize = 14;
colormap(lkcmap2)
Normalize_Clim(gcf,1)
savefig('Fig6c.fig')
close

%% Figure 6d: Electrode-based correlation with ERSP NAindex

figure
bandindex = find(strcmpi(fbands,'Alpha'));
nicecorrplot(nanmean(squeeze(restmeas.rel_bp.vals(bandindex,:,:))',1),nanmean(allmeas{bandindex}.naerspindex',1),{'Resting-state relative alpha power','Alpha nonadditivity'})
FixAxes(gca)
savefig('Fig6d.fig')
close

%% Figure 6e: Mediation model

opts = struct;
opts.display_mod = 1;
opts.display = 0;
opts.indvar = ['Rest' newline 'alpha'];
opts.depvar = ['Poststim' newline 'alpha'];
opts.mediator = ['Prestim' newline 'alpha'];
opts.sobelp = restmeas.rel_bp.mediation{4}.sobel.p;
mediationAnalysis0(double(nanmean(allmeas{bandindex}.erspindex(find(restmeas.rel_bp.prestim.stats{bandindex}.mask),:),1))',...
    double(squeeze(nanmean(restmeas.rel_bp.vals(bandindex,find(restmeas.rel_bp.prestim.stats{bandindex}.mask),:),2))),...
    double(nanmean(restmeas.prestimamp.rel{bandindex}(find(restmeas.rel_bp.prestim.stats{bandindex}.mask),:),1))',opts);
savefig('Fig6e.fig')
close
end
%% Combining the sub-figures

%% Figure 1

p = panel('no-manage-font')
p.pack('v',{50 25 25})
p(1).pack('h',{1/3 1/3 1/3})

p(1,1).pack('v',{50 50})
Fig1a = open('Fig1a.fig')
figaxes = findobj('Parent',Fig1a,'Type','axes');
p(1,1,1).select(figaxes(2))
FixAxes(p(1,1,1).axis,16)
p(1,1,2).select(figaxes(1))
FixAxes(p(1,1,2).axis,16)
close(Fig1a)

Fig1b = open('Fig1b.fig')
figaxes = findobj('Parent',Fig1b,'Type','axes');
p(1,2).select(figaxes)
FixAxes(p(1,2).axis,16)
close(Fig1b)

Fig1c = open('Fig1c.fig')
figaxes = findobj('Parent',Fig1c,'Type','axes');
p(1,3).select(figaxes)
FixAxes(p(1,3).axis,16)
close(Fig1c)

p(2).pack('h',repmat({1/settings.nfreqs},1,settings.nfreqs))

Fig1d = open('Fig1d.fig');
figaxes = findobj('Parent',Fig1d,'Type','Axes');

for c = 1:settings.nfreqs
    p(2,c).select(figaxes(settings.nfreqs+1-c));
    FixAxes(p(2,c).axis,12)
    set(gca,'TitleFontSizeMultiplier',1.1)
    if c ~= 1
        set(gca,'YLabel',[])
        
    end
end
close(Fig1d)

p(3).pack('h',[repmat({0.96/settings.nfreqs},1,settings.nfreqs) {0.03 0.01}])

Fig1e = open('Fig1e.fig');
figaxes = findobj('Parent',Fig1e,'Type','Axes');
cbar = findobj('Parent',Fig1e,'Type','ColorBar');


for c = 1:settings.nfreqs
    p(3,c).select(figaxes(settings.nfreqs+1-c));
    FixAxes(p(3,c).axis,12)
    set(gca,'TitleFontSizeMultiplier',1.1)
    if c ~= 1
        set(gca,'YLabel',[])
    end
end
p(3,settings.nfreqs+1).select(cbar)
colormap(lkcmap2)
close(Fig1e)

p.margin = [22 5 15 5];
p.de.margin = [5 5 5 5];
p(1).marginbottom = 22;
p(1,1,1).marginbottom = 18;
p(1,1).marginright = 18;
p(1,2).marginright = 18;
p(2).de.marginright = 10;
p(2).marginbottom = 20;

set(gcf,'Color','w')
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*2 pos(4)*2.5]) %make figure larger

AddFigureLabel(p(1,1,1).axis,'A')
AddFigureLabel(p(1,2).axis,'B')
AddFigureLabel(p(1,3).axis,'C')
AddFigureLabel(p(2,1).axis,'D')
AddFigureLabel(p(3,1).axis,'E')

savefig('Fig1.fig')
export_fig('Fig1.png','-m4')
save('Panel1.mat','p')
close

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
colormap(lkcmap2)
close(Fig2a)

p(2).pack()
p(2).pack({[0.02 0.02 0.3 0.3]})
Fig2b = open('Fig2b.fig')
figaxes = findobj('Parent',Fig2b,'Type','axes');
p(2,1).select()
nicecorrplot(nanmean(allmeas{plotband}.naerspindex,1),nanmean(allmeas{plotband}.ttverspindex,1),{'Pseudotrial-based nonadditivity','TTV-based nonadditivity'},'Plot','r')
FixAxes(p(2,1).axis,18)
p(2,2).select(figaxes(1))
set(p(2,2).axis,'FontSize',14)
colormap(lkcmap2)
close(Fig2b)

p.margin = [25 20 5 5];
p.de.margin = [5 5 5 5];
p(2).marginleft = 32;

AddFigureLabel(p(1,1).axis,'A')
AddFigureLabel(p(2,1).axis,'B')

set(gcf,'Color','w')
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*2 pos(4)*1.5]) %make figure larger
savefig('Fig2.fig')
export_fig('Fig2.png','-m4')
save('Panel2.mat','p')
close


%% Figure 3 already done
% figure
% p = panel('no-manage-font');
% p.pack('h',{32 32 32 2 2});
% Fig3 = open('Fig3.fig')
% figaxes = findobj('Parent',Fig3,'Type','axes')
% cbar = findobj('Parent',Fig3,'Type','Colorbar');
% p(1).select(figaxes(3))
% p(2).select(figaxes(2))
% p(3).select(figaxes(1))
% p(4).select(cbar)
%



%% Figure 4

figure
p = panel('no-manage-font');
p.pack('h',{50 50})
p(1).pack()
p(1).pack({[0.05 0.05 0.3 0.3]})
Fig4a = open('Fig4a.fig')
figaxes = findobj('Parent',Fig4a,'Type','axes')
p(1,1).select(figaxes(2));
FixAxes(p(1,1).axis,18)
p(1,2).select(figaxes(1));
set(p(1,2).axis,'FontSize',14)
close(Fig4a)
colormap(lkcmap2)

%FixAxes(p(1).axis,18)
Fig4b = open('Fig4b.fig')
figaxes = findobj('Parent',Fig4b,'Type','axes')
p(2).pack()
p(2).pack({[0.05 0.05 0.3 0.3]})
p(2,1).select();
nicecorrplot(nanmean(allmeas{1}.naerpindex,1)',nanmean(allmeas{1}.ttvindex,1)',{'Pseudotrial-based nonadditivity','TTV-based nonadditivity'},'Plot','r')
FixAxes(p(2,1).axis,18)
p(2,2).select(figaxes(1));
set(p(2,2).axis,'FontSize',14)
colormap(lkcmap2)
close(Fig4b)

p.marginbottom = 25;
p.marginleft = 30;
p(2).marginleft = 35;
p(1,2).margin = [5 5 5 5];
AddFigureLabel(p(1,1).axis,'A')
AddFigureLabel(p(2,1).axis,'B')



set(gcf,'Color','w')
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*2 pos(4)*1.5]) %make figure larger

savefig('Fig4.fig')
export_fig('Fig4.png','-m4')
save('Panel4.mat','p')
clear p
close

%% Figure 5

figure
p = panel('no-manage-font');
p.pack('v',{50 30 20});
p(1).pack(3,settings.nfreqs);
Fig5a = open('Fig5a.fig');
figaxes = findobj('Parent',Fig5a,'Type','axes');
figaxes = rot90(rot90(reshape(figaxes,3,settings.nfreqs)));
for c = 1:3
    for cc = 1:settings.nfreqs
        p(1,c,cc).select(figaxes(c,cc));
        FixAxes(p(1,c,cc).axis,12);
    end
end
close(Fig5a)

Fig5b = open('Fig5b_1.fig')
figaxes = findobj('Parent',Fig5b,'Type','axes');
p(2).select(figaxes);
FixAxes(p(2).axis,18);
close(Fig5b);

Fig5b = open('Fig5b_2.fig')
figaxes = findobj('Parent',Fig5b,'Type','axes');
figaxes = flipud(figaxes);
cbar = findobj('Parent',Fig5b','Type','ColorBar')
p(3).pack('h',{95 3 2});
p(3,1).pack('h',settings.nfreqs);
for c = 1:settings.nfreqs
    p(3,1,c).select(figaxes(c));
end
p(3,2).select(cbar)
colormap(lkcmap2)
close(Fig5b)

% Fig3c = open('Fig5c.fig')
% figaxes = findobj('Parent',Fig3c,'Type','axes');
% p(4).pack('h',{50 50})
% p(4,1).select(figaxes(2))
% FixAxes(p(4,1).axis,16)
% p(4,2).select(figaxes(1))
% FixAxes(p(4,2).axis,16)
p.margin = [28 8 8 10];
p.de.margin = [5 5 5 5];
p(1).margin = [5 18 5 5];
p(1).de.margin = [9 9 5 5];
p(2).marginbottom = 18;

AddFigureLabel(p(1,1,1).axis,'A');
AddFigureLabel(p(2).axis,'B');
AddFigureLabel(p(3,1,1).axis,'C');


set(gcf,'Color','w')
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*1.5 pos(4)*3]); %make figure larger
savefig('Fig5.fig')
export_fig('Fig5.png','-m4')
%p.export('Fig6','-rp')
save('Panel5.mat','p')
clear p
close

%% Figure 6
if isfield(settings,'rest')
figure
p = panel('no-manage-font');
p.pack('v',{26 26 48})
p(1).pack('h',{1/3 2/3})
p(1,2).pack(2,3)
Fig6a = open('Fig6a.fig');
figaxes = findobj('Parent',Fig6a,'Type','axes');
cbar = findobj('Parent',Fig6a,'Type','ColorBar');
%figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
plotindex = [(settings.nfreqs-1):-1:settings.nfreqs/2 ; (settings.nfreqs/2-1):-1:0];
for c = 1:2
    for cc = 1:settings.nfreqs/2
        if c == 2 && cc == settings.nfreqs/2
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
plotindex = [(settings.nfreqs-1):-1:settings.nfreqs/2 ; (settings.nfreqs/2-1):-1:0];
%figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
for c = 1:2
    for cc = 1:settings.nfreqs/2
        if c == 2 && cc == settings.nfreqs/2
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

set(gcf,'Color','w')
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*1.5 pos(4)*4]) %make figure larger
savefig('Fig6.fig')
export_fig('Fig6.png','-m4')
%p.export('Fig6','-rp')
save('Panel6.mat','p')
clear p
close
end
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

