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

%% Figure 0a. Nonadditivity of ERP - schematic

if strcmpi(settings.datatype,'MEG')
    plotsensor = find(strcmpi('A1',settings.datasetinfo.label));
else
    plotsensor = find(strcmpi('T7',settings.datasetinfo.label)); %just for Ivar data - where NA is strongest for schematic
end
bandindex = find(strcmpi('Broadband',fbands));
% t = linspace(-400,100,500);
% sig_low = sin(2*pi*t/100+rand*100).*[linspace(0.25,5,250) cos(2*pi*[0:249]/1000)];
% sig_high = sin(2*pi*t/100+rand*100).*[linspace(0.25,5,250) cos(2*pi*0:249/1000)];
% plot(sin(2*pi*t/100).*[linspace(0.25,5,250) cos(2*pi*0:249/1000)],'k','LineWidth',1.5);
subplot(3,1,1)
timeindex = [(settings.real.poststim(1)-settings.srate*0.4):(settings.real.poststim(1)+settings.srate*0.3)];
t = linspace(-400,300,length(timeindex));
plot(t,mean(allmeas{bandindex}.nadderp.raw.real(plotsensor,timeindex,1,:),4),'Color',[0.5 0.5 1],'LineWidth',0.5);
hold on
plot(t,mean(allmeas{bandindex}.nadderp.raw.real(plotsensor,timeindex,2,:),4),'Color',[1 0.5 0.5],'LineWidth',0.5);
xlabel('Time (ms)')
%legend({'Prestim low','Prestim High','Envelope of prestim low','Envelope of prestim high'});
title('Median split')
ylim = get(gca,'YLim');
line([0 0],ylim,'Color','k','LineWidth',2,'HandleVisibility','off')
patch([-length(settings.real.prestim)*1000/settings.srate 0 0 -length(settings.real.prestim)*1000/settings.srate],...
    [ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
set(gca,'YLim',ylim)
FixAxes(gca)
set(gca,'TitleFontSizeMultiplier',1.1)


subplot(3,1,2)
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
plot(t,mean(allmeas{1}.nadderp.real(plotsensor,:,1,:),4),'b','LineWidth',1.5)
plot(t,mean(allmeas{1}.nadderp.pseudo(plotsensor,:,1,:),4),'b--','LineWidth',1.5)
plot(t,mean(allmeas{1}.nadderp.pseudo(plotsensor,:,2,:),4),'r--','LineWidth',1.5)
plot(t,mean(allmeas{1}.nadderp.real(plotsensor,:,2,:),4),'r','LineWidth',1.5)
legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
%title('TTV')
FixAxes(gca)


subplot(3,1,3)
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
plot(t,mean(mean(allmeas{1}.nadderp.real(plotsensor,:,1,:),4),1)...
    -mean(mean(allmeas{1}.nadderp.pseudo(plotsensor,:,1,:),4),1),'b','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.nadderp.real(plotsensor,:,2,:),4),1)-...
    mean(mean(allmeas{1}.nadderp.pseudo(plotsensor,:,2,:),4),1),'r','LineWidth',1.5)
legend({'Corrected prestim low','Corrected prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
%title('Corrected TTV')
FixAxes(gca)

savefig('Fig0a.fig')


%% Figure 0b: Results of median split by band

%plotindex{1} = [3 4 5 8 9 10];
%plotindex{2} = [13 14 15 18 19 20];
figure
for c = 1:settings.nfreqs
    subplot(2,3,c)
    t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
    hold on
    plot(t,mean(mean(allmeas{c}.nadderp.real(:,:,1,:),4),1)-...
        mean(mean(allmeas{c}.nadderp.pseudo(:,:,1,:),4),1),'b','LineWidth',1.5)
    
    plot(t,mean(mean(allmeas{c}.nadderp.real(:,:,2,:),4),1)-...
        mean(mean(allmeas{c}.nadderp.pseudo(:,:,2,:),4),1),'r','LineWidth',1.5)
    title(fbands{c})
    %legend({'Corrected prestim low','Corrected prestim high'})
    xlabel('Time (s)')
    ylabelunits(settings)
    FixAxes(gca)
    set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
end
savefig(['Fig_0b.fig'])
close

%% Figure 0c. Significance of NAindex for ERP

figure
for c = 1:settings.nfreqs
    ax(c) = subplot(4,3,c);
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,median(allmeas{c}.naerpindex,2).*...
            (alloutputs.erp.stats{c}.mask),settings.datasetinfo.label);
    else
        topoplot(median(allmeas{c}.naerpindex,2).*...
            (alloutputs.erp.stats{c}.mask),settings.layout);
    end
    title(fbands{c})
    set(gca,'TitleFontSizeMultiplier',1.3)
end
savefig(['Fig0c.fig'])
close

%% Figure 1a. Raw TTV, Amplitude, ITC time courses

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
    
    ax(c+settings.nfreqs) = subplot(3,settings.nfreqs,c+settings.nfreqs);
    stdshade(t,squeeze(mean(allmeas{c}.raw.ersp(:,poststim_real,:),1)),'r',0.15,1)
    hold on
    stdshade(t,squeeze(mean(allmeas{c}.raw.ersp(:,poststim_pseudo,:),1)),'r--',0.15,1)
    if c == 1
        ylabel('Hilbert amplitude')
    end
    FixAxes(gca)
    
    ax(c+2*settings.nfreqs) = subplot(3,settings.nfreqs,c+2*settings.nfreqs);
    stdshade(t,squeeze(mean(allmeas{c}.raw.itc(:,poststim_real,:),1)),'g',0.15,1)
    hold on
    stdshade(t,squeeze(mean(allmeas{c}.raw.itc(:,poststim_pseudo,:),1)),'g--',0.15,1)
    xlabel('Time (s)')
    
    if c == 1
        ylabel('Inter-trial coherence')
    end
    FixAxes(gca)
    
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

savefig('Fig1a.fig')
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
    FixAxes(ax(c))
end

savefig('Fig1b.fig')
close

%% Figure 1c. Time course similarity topoplots

figure

for c = 1:6
    subplot(2,3,c)
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,(mean(alloutputs.itc.distrealreal(:,:,c),2)-mean(alloutputs.ersp.distrealreal(:,:,c),2)).*...
            (alloutputs.ttv.diststats{c}.mask > 0),settings.datasetinfo.label);
    else
        topoplot((mean(alloutputs.itc.distrealreal(:,:,c),2)-mean(alloutputs.ersp.distrealreal(:,:,c),2)).*...
            (alloutputs.ttv.diststats{c}.mask > 0),settings.layout);
    end
    title(fbands{c})
    set(gca,'TitleFontSizeMultiplier',1.3)
end
cbar = colorbar('peer',gca,'Position',...
    [0.90287 0.37857 0.02429 0.3071429],'FontSize',12);
cbar.Label.String = 'Euclidean Distance Difference';
cbar.Label.FontSize = 14;
savefig('Fig1c.fig')
close


%% Figure 2a. Schematic of correlation between TTVindex and ERSPindex, ITCindex

%designed for 6 frequency bands

figure
subplot(4,4,[5 6 9 10])
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
set(gca,'FontSize',14,'position',p+[-0.02 0 -0.04 0],'TitleFontSizeMultiplier',1.3,...
    'XLim',[0 0.8])
FixAxes(gca)


subplot(4,4,[3 4 7 8])
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
set(gca,'FontSize',14,'position',p+[0.05 0 0 0],'TitleFontSizeMultiplier',1.3)
FixAxes(gca)

subplot(4,4,[11 12 15 16])
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
set(gca,'FontSize',14,'position',p+[0.05 0 0 0],'TitleFontSizeMultiplier',1.3)
FixAxes(gca)

% Create arrow
annotation(gcf,'arrow',[0.441908713692946 0.518950437317784],...
    [0.553376657824934 0.752711496746204]);

% Create arrow
annotation(gcf,'arrow',[0.444606413994169 0.510204081632653],...
    [0.458869848156182 0.253796095444686]);

% Create textbox
annotation(gcf,'textbox',...
    [0.37227 0.64595 0.07 0.07],...
    'String',{'TTV index'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.885375 0.87194 0.1 0.07],...
    'String','ERSP index',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.89553125 0.4368455 0.1 0.07],...
    'String','ITC index',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');

savefig('Fig2a.fig')
close

%% Figure 2b. Scatter plot of index correlation for ERSP, ITC broadband

figure
subplot(1,2,1)
nicecorrplot(mean(allmeas{1}.ttvindex(:,:),1),mean(allmeas{1}.erspindex(:,:),1),{'TTV Index','ERSP Index'});
FixAxes(gca)
subplot(1,2,2)
nicecorrplot(mean(allmeas{1}.ttvindex(:,:),1),mean(allmeas{1}.itcindex(:,:),1),{'TTV Index','ITC Index'});
FixAxes(gca)
savefig('Fig2b.fig')

%% Figure 2c. Topoplots of index correlation
figure
plotindex1 = [1 2 3 7 8 9];
plotindex2 = [4 5 6 10 11 12];
for c = 1:settings.nfreqs
    subplot(2,6,plotindex1(c))
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,alloutputs.ersp.indexr(:,c).*(alloutputs.ersp.corrstats{c}.mask > 0),settings.datasetinfo.label);
    else
        topoplot(alloutputs.ersp.indexr(:,c).*(alloutputs.ersp.corrstats{c}.mask > 0),settings.layout);
    end
    title(fbands{c})
    
    subplot(2,6,plotindex2(c))
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

savefig('Fig2c.fig')
close

% Create textbox
% annotation(gcf,'textbox',[0.065 0.67 0.015 0.045],'String','a',...
%     'FontSize',24,'FitBoxToText','off','EdgeColor','none');
%
% % Create textbox
% annotation(gcf,'textbox',[0.31 0.8625 0.0180 0.06285],'String',{'b'},...
%     'FontSize',24,'FitBoxToText','off','EdgeColor','none');
%
% % Create textbox
% annotation(gcf,'textbox',[0.56 0.8625 0.018 0.06285],'String','c',...
%     'FontSize',24,'FitBoxToText','off','EdgeColor','none');
%
% % Create textbox
% annotation(gcf,'textbox',[0.56 0.445 0.0180 0.06285],'String','d',...
%     'FontSize',24,'FitBoxToText','off','EdgeColor','none');



%% Figure 3a. Split schematic, real vs pseudo, and corrected real trials

if strcmpi(settings.datatype,'MEG')
    plotsensor = find(strcmpi('A1',settings.datasetinfo.label));
else
    plotsensor = find(strcmpi('CPz',settings.datasetinfo.label));
end
bandindex = find(strcmpi('Alpha',fbands));
% t = linspace(-400,100,500);
% sig_low = sin(2*pi*t/100+rand*100).*[linspace(0.25,5,250) cos(2*pi*[0:249]/1000)];
% sig_high = sin(2*pi*t/100+rand*100).*[linspace(0.25,5,250) cos(2*pi*0:249/1000)];
% plot(sin(2*pi*t/100).*[linspace(0.25,5,250) cos(2*pi*0:249/1000)],'k','LineWidth',1.5);
subplot(3,2,1)
timeindex = [(settings.real.poststim(1)-settings.srate*0.4):(settings.real.poststim(1)+settings.srate*0.3)];
t = linspace(-400,300,length(timeindex));
plot(t,mean(allmeas{bandindex}.naddersp.amp.raw.real(plotsensor,timeindex,1,:),4),'Color',[0.5 0.5 1],'LineWidth',0.5);
hold on
plot(t,mean(allmeas{bandindex}.naddersp.amp.raw.real(plotsensor,timeindex,2,:),4),'Color',[1 0.5 0.5],'LineWidth',0.5);
plot(t,abs(hilbert(mean(allmeas{bandindex}.naddersp.amp.raw.real(plotsensor,timeindex,1,:),4))),'b','LineWidth',1.5);
plot(t,abs(hilbert(mean(allmeas{bandindex}.naddersp.amp.raw.real(plotsensor,timeindex,2,:),4))),'r','LineWidth',1.5);
xlabel('Time (ms)')
ylabel('Voltage (uv)')
%legend({'Prestim low','Prestim High','Envelope of prestim low','Envelope of prestim high'});
title('Envelope-based split')
ylim = get(gca,'YLim');
line([0 0],ylim,'Color','k','LineWidth',2,'HandleVisibility','off')
patch([-200 0 0 -200],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
set(gca,'YLim',ylim)
FixAxes(gca)
set(gca,'TitleFontSizeMultiplier',1.1)


subplot(3,2,2)
timeindex = [(settings.real.poststim(1)-settings.srate*0.4):(settings.real.poststim(1)+settings.srate*0.3)];
t = linspace(-400,300,length(timeindex));
plot(t,mean(mean(real(allmeas{bandindex}.naddersp.cos.raw.real(plotsensor,timeindex,1,:)),4),1),'Color',[0.5 0.5 1],'LineWidth',0.5);
hold on
plot(t,mean(mean(real(allmeas{bandindex}.naddersp.cos.raw.real(plotsensor,timeindex,2,:)),4),1),'Color',[1 0.5 0.5],'LineWidth',0.5);
%plot(t,cos(angle(hilbert(mean(allmeas{bandindex}.naddersp.cos.raw.real(plotsensor,timeindex,1,:),4)))),'b','LineWidth',1.5);
%plot(t,cos(angle(hilbert(mean(allmeas{bandindex}.naddersp.cos.raw.real(plotsensor,timeindex,2,:),4)))),'r','LineWidth',1.5);
%legend({'Prestim low','Prestim High','Cos of phase of prestim low','Cos of phase of prestim high'})
xlabel('Time (ms)')
ylabel('Voltage (uv)')
title('Phase-based split')
ylim = get(gca,'YLim');
line([0 0],ylim,'Color','k','LineWidth',2)
%patch([-200 0 0 -200],[ylim(1) ylim(1) ylim(2) ylim(2)],[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');
set(gca,'YLim',ylim)
FixAxes(gca)
set(gca,'TitleFontSizeMultiplier',1.1)


subplot(3,2,[3 4])
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,1,:),4),1),'b','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,1,:),4),1),'b--','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,2,:),4),1),'r--','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,2,:),4),1),'r','LineWidth',1.5)
legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
%title('TTV')
FixAxes(gca)


subplot(3,2,[5 6])
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,1,:),4),1)...
    -mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,1,:),4),1),'b','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,2,:),4),1)-...
    mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,2,:),4),1),'r','LineWidth',1.5)
legend({'Corrected prestim low','Corrected prestim high'})
xlabel('Time (s)')
ylabelunits(settings)
%title('Corrected TTV')
FixAxes(gca)

savefig('Fig3a.fig')


%% Figures 3b and 3c: Results for power- and phase-based split of TTV

splitmethod = {'amp','cos'};
%plotindex{1} = [3 4 5 8 9 10];
%plotindex{2} = [13 14 15 18 19 20];
for q = 1:2
    figure
    for c = 1:settings.nfreqs
        subplot(1,6,c)
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
        FixAxes(gca)
        set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
    end
    savefig(['Fig_3bcde_' num2str(q) '.fig'])
    close
end

%% Figures 3d and 3e: Results for power- and phase-based split of ERSP

splitmethod = {'amp','cos'};
% plotindex{1} = [3 4 5 8 9 10];
% plotindex{2} = [13 14 15 18 19 20];
for q = 1:2
    figure
    for c = 1:settings.nfreqs
        subplot(1,6,c)
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
        FixAxes(gca)
        set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
    end
    savefig(['Fig_3bcde_' num2str(q+2) '.fig'])
    close
end

%% Figure 4a. Schematic of nonadditivity index

%subplot(3,2,[5 6])
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
x2 = [t(settings.aucindex), fliplr(t(settings.aucindex))];
inBetween = [(mean(mean(allmeas{1}.naddttv.amp.real(:,:,1,:),4),1)...
    -mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,1,:),4),1)), ...
    fliplr(mean(mean(allmeas{1}.naddttv.amp.real(:,:,2,:),4),1)-...
    mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,2,:),4),1))];
fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');

hold on
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,1,:),4),1)...
    -mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,1,:),4),1),'b','LineWidth',1.5)
plot(t,mean(mean(allmeas{1}.naddttv.amp.real(:,:,2,:),4),1)-...
    mean(mean(allmeas{1}.naddttv.amp.pseudo(:,:,2,:),4),1),'r','LineWidth',1.5)
legend({'Corrected prestim low','Corrected prestim high','Nonadditivity Index'})
xlabel('Time (s)')
ylabelunits(settings)
title('Nonadditivity index')
FixAxes(gca)
savefig('Fig4a.fig')
close


%% Figure 4bcde. Significance of NAindex for ERSP and TTV
bcde = {'b','c','d','e'};
for q = 1:2
    figure
    for c = 1:settings.nfreqs
        ax(c) = subplot(1,6,c);
        if strcmpi(settings.datatype,'MEG')
            ft_topoplot_vec(settings.layout,median(allmeas{c}.nattvindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ttv.stats{c}.mask),settings.datasetinfo.label);
        else
            topoplot(median(allmeas{c}.nattvindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ttv.stats{c}.mask),settings.layout);
        end
        title(fbands{c})
        set(gca,'TitleFontSizeMultiplier',1.3)
    end
    Normalize_Clim(gcf)
    
    %Aesthetics
    cbar = colorbar('peer',ax(1),'FontSize',12);
    cbar.Label.String = 'Nonadditivity Index';
    cbar.Label.FontSize = 14;
    
    savefig(['Fig4bcde_' num2str(q)])
    close
end

for q = 1:2
    for c = 1:settings.nfreqs
        ax(c) = subplot(1,6,c);
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
    cbar = colorbar('peer',ax(1),'FontSize',12);
    cbar.Label.String = 'Nonadditivity Index';
    cbar.Label.FontSize = 14;
    
    %     % Create textbox
    %     annotation(gcf,'textbox',...
    %         [0.16178 0.861554 0.01947 0.07283],...
    %         'String',bcde{(q-1)*2+1},...
    %         'FontSize',24,...
    %         'FitBoxToText','off',...
    %         'EdgeColor','none');
    %
    %     % Create textbox
    %     annotation(gcf,'textbox',...
    %         [0.16178 0.423188 0.01947 0.07283],...
    %         'String',bcde{2*q},...
    %         'FontSize',24,...
    %         'FitBoxToText','off',...
    %         'EdgeColor','none');
    
    savefig(['Fig4bcde_' num2str(q+2)])
    close
end


%% Figure 5a: Topoplots of resting state correlations with ERSPindex

figure

for c = 2:6
    ax(c) = subplot(2,3,c-1);
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,restmeas.rel_bp.index.r.subject(:,c).*...
            (restmeas.rel_bp.index.stats{c}.mask),settings.datasetinfo.label);
    else
        topoplot(restmeas.rel_bp.index.r.subject(:,c).*...
            (restmeas.rel_bp.index.stats{c}.mask),settings.layout);
    end
    title(settings.tfparams.fbandnames{c})
end
set(gca,'TitleFontSizeMultiplier',1.3)
cbar = colorbar('peer',gca,'FontSize',12);
cbar.Label.String = "Spearman's rho";
cbar.Label.FontSize = 14;
savefig('Fig5a.fig')

%% Figure 5b: Electrode-based correlation with ERSPindex
figure
bandindex = find(strcmpi(fbands,'Alpha'));
nicecorrplot(mean(squeeze(restmeas.rel_bp.vals(bandindex,:,:)),2),mean(allmeas{bandindex}.erspindex,2),{'Resting-state relative alpha power','Alpha ERSPindex'})
FixAxes(gca)
savefig('Fig5b.fig')

%% Figure 5c: Topoplots of resting state correlations with ERSP NAindex
for c = 2:6
    ax(c) = subplot(2,3,c-1);
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,restmeas.rel_bp.naindex.r.subject(:,c).*...
            (restmeas.rel_bp.naindex.stats{c}.mask),settings.datasetinfo.label);
    else
        topoplot(restmeas.rel_bp.naindex.r.subject(:,c).*...
            (restmeas.rel_bp.naindex.stats{c}.mask),settings.layout);
    end
    title(settings.tfparams.fbandnames{c})
end
set(gca,'TitleFontSizeMultiplier',1.3)
cbar = colorbar('peer',gca','FontSize',12);
cbar.Label.String = "Spearman's rho";
cbar.Label.FontSize = 14;
savefig('Fig5c.fig')

%% Figure 5d: Electrode-based correlation with ERSP NAindex

figure
bandindex = find(strcmpi(fbands,'Alpha'));
nicecorrplot(mean(squeeze(restmeas.rel_bp.vals(bandindex,:,:)),2),mean(allmeas{bandindex}.naerspindex.amp,2),{'Resting-state relative alpha power','Alpha Nonadditivity Index'})
FixAxes(gca)
savefig('Fig5d.fig')

%% Figure 5e: Mediation model
figure
mediationAnalysis0(double(mean(allmeas{bandindex}.erspindex(find(restmeas.rel_bp.prestim.stats{bandindex}.mask),:),1))',...
    double(squeeze(mean(restmeas.rel_bp.vals(bandindex,find(restmeas.rel_bp.prestim.stats{bandindex}.mask),:),2))),...
    double(mean(restmeas.prestimamp.rel{bandindex}(find(restmeas.rel_bp.prestim.stats{bandindex}.mask),:),1))');
savefig('Fig5e.fig')

%
% figure
% if isfield(settings,'rest')
%     subplot(2,5,[1 2 6 7])
%     nicecorrplot(squeeze(mean(restmeas.rel_bp.vals(4,:,:),2)),squeeze(mean(allmeas{4}.erspindex,1)),{'Relative resting-state band-limited power','Poststimulus ERSPindex'});
%     FixAxes(gca)
%     plotindex = [3 4 5 8 9];
%     for c = 2:settings.nfreqs
%         ax(c) = subplot(2,5,plotindex(c-1));
%         if strcmpi(settings.datatype,'MEG')
%             ft_topoplot_vec(settings.layout,restmeas.rel_bp.index.r(:,c).*...
%                 (restmeas.rel_bp.index.stats{c}.mask),settings.datasetinfo.label);
%         else
%             topoplot(restmeas.rel_bp.index.r(:,c).*...
%                 (restmeas.rel_bp.index.stats{c}.mask),settings.layout);
%         end
%         title(settings.tfparams.fbandnames{c})
%         set(gca,'TitleFontSizeMultiplier',1.3)
%     end
%
%     pos = get(gcf,'position');
%     set(gcf,'position',[pos(1:2)/2 pos(3)*2 pos(4)*1.5])
%     % Create colorbar
%     cbar = colorbar('peer',ax(2),'Position',...
%         [0.801837 0.148936 0.020385668 0.268461],'FontSize',12);
%     cbar.Label.String = 'Spearman Correlation'
%     cbar.Label.FontSize = 14;
%
%     % Create textbox
%     annotation(gcf,'textbox',...
%         [0.0520725 0.870 0.01824 0.05816],...
%         'String','a',...
%         'FontSize',24,...
%         'FitBoxToText','on',...
%         'EdgeColor','none');
%
%     % Create textbox
%     annotation(gcf,'textbox',...
%         [0.43410 0.870 0.01824 0.0581],...
%         'String','b',...
%         'FontSize',24,...
%         'FitBoxToText','on',...
%         'EdgeColor','none');
%     savefig(['Fig_5a.fig'])

%% Combining the sub-figures

%% Figure 0

figure
p = panel();
p.pack('h',{1/3 2/3})
p(2).pack(2,1)

p(1).pack(3,1)
Fig0a = open('Fig0a.fig')
figaxes = findobj('Parent',Fig1a,'Type','axes');
p(1,1).select(figaxes(3))
FixAxes(p(1,1).axis)
p(1,2).select(figaxes(2))
FixAxes(p(1,2).axis)
p(1,3).select(figaxes(1))
FixAxes(p(1,3).axis)
close(Fig0a)

p(2,1).pack(2,3)
Fig0b = open('Fig0b.fig');
figaxes = findobj('Parent',Fig0b,'Type','axes');
figaxes = rot90(rot90(reshape(figaxes,2,3)));
for c = 1:2
    for cc = 1:3
        p(2,1,c,cc).select(figaxes(c,cc));
        FixAxes(p(2,1,c,cc).axis);
    end
end
close(Fig0b)

p(2,2).pack('h',{85 15}) %leaving room for colorbar
Fig0c = open('Fig0c.fig');
figaxes = findobj('Parent',Fig0c,'Type','axes');
figaxes = rot90(rot90(reshape(figaxes,2,3)));
for c = 1:2
    for cc = 1:3
        p(2,2,1,c,cc).select(figaxes(c,cc));
        FixAxes(p(2,2,1,c,cc).axis);
    end
end
NormalizeClim(gcf)
cbar = colorbar('peer',p(2,2,1,2,3).axis,'FontSize',12);
cbar.Label.String = 'Nonadditivity index';
cbar.Label.FontSize = 14;
p(2,2,2).select(cbar)
close(Fig0c)

savefig('Fig0.fig')
p.export('Fig0','-rp')
save('Panel0.mat','p')
clear p
close



%% Figure 1

figure
p = panel('no-manage-font');
p.pack('h',{60 40})
p(2).pack('v',{50 50})

p(1).pack(3,6)
Fig1a = open('Fig1a.fig');
figaxes = findobj('Parent',Fig1a,'Type','axes');
figaxes = rot90(rot90(reshape(figaxes,3,6)));
for c = 1:3
    for cc = 1:6
        p(1,c,cc).select(figaxes(c,cc));
        FixAxes(p(1,c,cc).axis);
    end
end
%p.marginleft = 30;
AddFigureLabel(p(1,1,1).position,'a','panel')
close(Fig1a)

p(2,1).pack(2,3)
Fig1b = open('Fig1b.fig')
figaxes = findobj('Parent',Fig1b,'Type','axes');
figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
for c = 1:2
    for cc = 1:3
        p(2,1,c,cc).select(figaxes(c,cc));
        FixAxes(p(2,1,c,cc).axis);
    end
end
AddFigureLabel(p(2,1,1,1).position,'b','panel')
close(Fig1b)

p(2,2).pack('h',{93 7});
p(2,2,1).pack(2,3)
Fig1c = open('Fig1c.fig')
figaxes = findobj('Parent',Fig1c,'Type','axes');
cbar = findobj('Parent',Fig1c,'Type','ColorBar');
figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
for c = 1:2
    for cc = 1:3
        p(2,2,1,c,cc).select(figaxes(c,cc));
        FixAxes(p(2,2,1,c,cc).axis);
    end
end
AddFigureLabel(p(2,2,1,1,1).position,'c','panel')
p(2,2,2).select(cbar)
%p(2,2,1,2).margin = [20 100 50 50];
close(Fig1c)

%Fix margins
p.margin = [30 30 30 15];
p.de.margin = [15 15 5 5];
p(1).marginright = 30;

savefig('Fig1.fig')
p.export('Fig1','-rp')
save('Panel1.mat','p')
clear p
close

%% Figure 2
figure

p = panel('no-manage-font');
p.margin = [30 30 30 15];
p.de.margin = [15 15 5 5];
p.pack('v',{2/3 []})
p(1).pack('h',{2/3 []})
p(1,1).pack('h',{50 50})
p(1,1,1).pack('v',{20 60 20})
p(1,1,2).pack('v',{50 50})

Fig2a = open('Fig2a.fig')
figaxes = findobj('Parent',Fig2a,'Type','axes');
p(1,1,1,2).select(figaxes(3))
FixAxes(p(1,1,1,2).axis)
p(1,1,2,1).select(figaxes(2))
FixAxes(p(1,1,2,1).axis)
p(1,1,2,2).select(figaxes(1))
FixAxes(p(1,1,2,2).axis)
AddFigureLabel(p(1,1,1,1).position,'a','panel')
close(Fig2a)

p(1,2).pack('v',{50 50})
Fig2b = open('Fig2b.fig')
figaxes = findobj('Parent',Fig2b,'Type','axes');
p(1,2,1).select(figaxes(2));
FixAxes(p(1,2,1).axis)
p(1,2,2).select(figaxes(1));
FixAxes(p(1,2,2).axis)
AddFigureLabel(p(1,2,1).position,'b','panel')
close(Fig2b)

p(2).pack('h',{97 3})
p(2,1).pack(2,6)
Fig2c = open('Fig2c.fig');
figaxes = findobj('Parent',Fig2c,'Type','axes');
cbar = findobj('Parent',Fig2c,'Type','ColorBar');
figaxes = fliplr(flipud(transpose(reshape(figaxes,6,2))));
for c = 1:2
    for cc = 1:6
        p(2,1,c,cc).select(figaxes(c,cc))
        FixAxes(p(2,1,c,cc).axis)
    end
end
AddFigureLabel(p(2,1,1,1).position,'c','panel')
p(2,2).select(cbar)
close(Fig2c)

%Fix margins
p.margin = [30 30 30 15];
p.de.margin = [15 15 5 5];
p(1).marginbottom = 30;
p(1,1).marginright = 30;
p(1,1,1).marginright = 25;

savefig('Fig2.fig')
p.export('Fig2','-rp')
save('Panel2.mat','p')
clear p
close

%% Figure 3
p = panel('no-manage-font');
p.pack('h',{30 70})
p(1).pack('v',{1/3 1/3 1/3})
p(1,1).pack('h',{50 50})
Fig3a = open('Fig3a.fig');
figaxes = findobj('Parent',Fig3a,'Type','axes');
p(1,1,1).select(figaxes(4))
FixAxes(p(1,1,1).axis)
p(1,1,2).select(figaxes(3))
FixAxes(p(1,1,2).axis)
p(1,2).select(figaxes(2))
FixAxes(p(1,2).axis)
p(1,3).select(figaxes(1))
FixAxes(p(1,3).axis)
AddFigureLabel(p(1,1,1).position,'a')
close(Fig3a)

p(2).pack('v',{25 25 25 25})
%layout = {[1 1], [1 2], [2 1], [2 2]};
for c = 1:4
    p(2,c).pack(1,6)
    Fig3bcde = open(['Fig_3bcde_' num2str(c) '.fig'])
    figaxes = findobj('Parent',Fig3bcde,'Type','axes');
    figaxes = flipud(figaxes);
    for cc = 1:6
        p(2,c,1,cc).select(figaxes(cc))
        FixAxes(p(2,c,1,cc).axis)
        if c > 1
            set(gca,'Title',[])
        end
        if cc > 1
            set(gca,'YLabel',[])
        end
        if c < 4
            set(gca,'XLabel',[])
        end
    end
    AddFigureLabel(p(2,c,1,1).position,bcde{c},'panel')
    close(Fig3bcde)
end

%Fix margins
p.margin = [20 20 10 10];
p.de.margin = [10 10 5 5];
p(1).marginright = 25;
p(1).de.marginbottom = 20;

savefig('Fig3.fig')
p.export('Fig3','-rp')
save('Panel3.mat','p')
clear p
close

%% Figure 4

figure
p = panel('no-manage-font');
p.pack('v',{25 75})
p(1).pack('h',{25 50 25})
Fig4a = open('Fig4a.fig');
figaxes = findobj('Parent',Fig4a,'Type','axes');
p(1,2).select(figaxes)
AddFigureLabel(p(1,2).position,'a')
FixAxes(p(1,2).axis)
close(Fig4a)

p(2).pack('h',{98 2})
p(2,1).pack('v',{25 25 25 25})
p(2,2).pack('v',{25 25 25 25})
for c = 1:4
    p(2,1,c).pack(1,6)
    Fig4bcde = open(['Fig4bcde_' num2str(c) '.fig'])
    figaxes = findobj('Parent',Fig4bcde,'Type','axes');
    cbar = findobj('Parent',Fig4bcde,'Type','ColorBar');
    figaxes = flipud(figaxes);
    for cc = 1:6
        p(2,1,c,1,cc).select(figaxes(cc))
        FixAxes(p(2,1,c,1,cc).axis)
        if c > 1
            set(gca,'Title',[])
        end
    end
    AddFigureLabel(p(2,1,1,1,1).position,bcde{c})
    p(2,2,c).select(cbar)
    close(Fig4bcde)
end

%Fix margins
p.margin = [15 15 22 10];
p.de.margin = [10 10 5 5];
p(1).marginbottom = 30;
p(2,1).marginleft = 15;

savefig('Fig4.fig')
p.export('Fig4','-rp')
save('Panel4.mat','p')
clear p
close

%% Figure 5

figure
p = panel('no-manage-font');
p.pack('v',{35 35 30})
p(1).pack('h',{2/3 1/3})
p(1,1).pack(2,6)
Fig5a = open('Fig5a.fig');
figaxes = findobj('Parent',Fig5a,'Type','axes');
cbar = findobj('Parent',Fig5a,'Type','ColorBar');
%figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
plotindex = [5 4 3; 2 1 0];
for c = 1:2
    for cc = 1:3
        if c == 2 && cc == 3
            p(1,1,c,cc).pack('h',{90 10});
            p(1,1,c,cc,2).select(cbar);
        else
        p(1,1,c,cc).select(figaxes(plotindex(c,cc)));
        FixAxes(p(1,1,c,cc).axis)
        end
    end
end
close(Fig5a)

Fig5b = open('Fig5b.fig');
figaxes = findobj('Parent',Fig5b,'Type','axes');
p(1,2).select(figaxes)
FixAxes(p(1,2).axis)
close(Fig5b)

p(2).pack('h',{2/3 1/3})
p(2,1).pack(2,6)
Fig5c = open('Fig5c.fig');
figaxes = findobj('Parent',Fig5c,'Type','axes');
cbar = findobj('Parent',Fig5c,'Type','ColorBar');
plotindex = [5 4 3; 2 1 0];
%figaxes = fliplr(flipud(transpose(reshape(figaxes,3,2))));
for c = 1:2
    for cc = 1:3
        if c == 2 && cc == 3
            p(2,1,c,cc).pack('h',{90 10});
            p(2,1,c,cc,2).select(cbar);
        else
        p(2,1,c,cc).select(figaxes(plotindex(c,cc)));
        FixAxes(p(2,1,c,cc).axis)
        end
    end
end
close(Fig5c)

Fig5d = open('Fig5d.fig');
figaxes = findobj('Parent',Fig5d,'Type','axes');
p(2,2).select(figaxes)
FixAxes(p(2,2).axis)
close(Fig5d)

Fig5e = open('Fig5e.fig');
figaxes = findobj('Parent',Fig5e,'Type','axes');
p(3).select(figaxes)
close(Fig5e)

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

