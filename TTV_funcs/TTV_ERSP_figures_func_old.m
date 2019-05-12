function TTV_ERSP_figures_func(settings)

fbands = settings.tfparams.fbandnames;

load([settings.outputdir '/' settings.datasetname '_calc.mat'])
load([settings.outputdir '/' settings.datasetname '_results.mat'])

if strcmpi(settings.datatype,'EEG')
    eeglab
end

cd([settings.outputdir '/' settings.datasetname '_figures'])

prestim_pseudo = settings.pseudo.prestim;
prestim_real = settings.real.prestim;
poststim_pseudo = settings.pseudo.poststim;
poststim_real = settings.real.poststim;

%% Figure 1. Raw TTV, Amplitude, ITC time courses

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
    ylabel('Amplitude')
    
    subplot(3,settings.nfreqs,c+2*settings.nfreqs)
    plot(t,mean(mean(allmeas{c}.raw.itc(:,poststim_real,:),3),1),'g','LineWidth',1.5)
    hold on
    plot(t,mean(mean(allmeas{c}.raw.itc(:,poststim_real,:),3),1),'g--','LineWidth',1.5)
    
    xlabel('Time (s)')
    ylabel('Inter-trial coherence')
    
end
savefig('Fig1.fig')
close

%% Figure 2a. Average TTV, ERSP, ITC time courses for each band
figure
for c = 1:settings.nfreqs
    subplot(2,4,c)
    t = linspace(0,max(settings.real.poststim)/settings.srate,length(settings.real.poststim));
    plot(t,mean(mean(allmeas{c}.ttv.real,3),1),'b','LineWidth',1.5)
    hold on
    plot(t,mean(mean(allmeas{c}.ersp.real,3),1),'r','LineWidth',1.5)
    plot(t,mean(mean(allmeas{c}.itc.real,3),1),'g','LineWidth',1.5)
    
    title(fbands{c})
    legend({'TTV','ERSP','ITC'})
    xlabel('Time (s)')
    ylabelunits(settings)
end

savefig('Fig2a.fig')
close

%% Figure 2b. Correlation between TTVindex and ERSPindex, ITCindex

%designed for 6 frequency bands

figure
subplot(7,4,[8 9 15 16])
t = linspace(0,max(settings.real.poststim)/settings.srate,length(settings.real.poststim));
plot(t,mean(mean(allmeas{1}.ttv.real,3),1),'b','LineWidth',1.5); %assumes 1 is broadband
hold on
plot(t,mean(mean(allmeas{1}.ttv.pseudo,3),1),'b--','LineWidth',1.5)
x2 = [t(settings.aucindex), fliplr(t(settings.aucindex))];
inBetween = [mean(mean(allmeas{1}.ttv.real(:,settings.aucindex,:),3),1), ...
    fliplr(mean(mean(allmeas{1}.ttv.pseudo(:,settings.aucindex,:),3),1))];
fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.3);
xlabel('Time (s)')
ylabelunits(settings)
    
subplot(7,4,[3 4 10 11])
t = linspace(0,max(settings.real.poststim)/settings.srate,length(settings.real.poststim));
plot(t,mean(mean(allmeas{1}.ersp.real,3),1),'b','LineWidth',1.5); %assumes 1 is broadband
hold on
plot(t,mean(mean(allmeas{1}.ersp.pseudo,3),1),'b--','LineWidth',1.5)
x2 = [t(settings.aucindex), fliplr(t(settings.aucindex))];
inBetween = [mean(mean(allmeas{1}.ersp.real(:,settings.aucindex,:),3),1), ...
    fliplr(mean(mean(allmeas{1}.ersp.pseudo(:,settings.aucindex,:),3),1))];
fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.3);
xlabel('Time (s)')
ylabelunits(settings)

subplot(7,4,[17 18 24 25])
t = linspace(0,max(settings.real.poststim)/settings.srate,length(settings.real.poststim));
plot(t,mean(mean(allmeas{1}.itc.real,3),1),'b','LineWidth',1.5); %assumes 1 is broadband
hold on
plot(t,mean(mean(allmeas{1}.itc.pseudo,3),1),'b--','LineWidth',1.5)
x2 = [t(settings.aucindex), fliplr(t(settings.aucindex))];
inBetween = [mean(mean(allmeas{1}.itc.real(:,settings.aucindex,:),3),1), ...
    fliplr(mean(mean(allmeas{1}.itc.pseudo(:,settings.aucindex,:),3),1))];
fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.3);
xlabel('Time (s)')
ylabelunits(settings)

plotindex1 = [5 6 7 12 13 14]; 
plotindex2 = [19 20 21 26 27 28];
for c = 1:settings.nfreqs
    subplot(7,4,plotindex1(c))
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,alloutputs.ersp.indexr(:,c).*(alloutputs.ersp.corrstats{c}.mask > 0),settings.datasetinfo.label);
    else
        topoplot(alloutputs.ersp.indexr(:,c).*(alloutputs.ersp.corrstats{c}.mask > 0),settings.layout);
    end
    title(fbands{c})
    
    subplot(4,4,plotindex2(c))
    if strcmpi(settings.datatype,'MEG')
        ft_topoplot_vec(settings.layout,alloutputs.itc.indexr(:,c).*(alloutputs.itc.corrstats{c}.mask > 0),settings.datasetinfo.label);
    else
        topoplot(alloutputs.itc.indexr(:,c).*(alloutputs.itc.corrstats{c}.mask > 0),settings.layout);
    end
    title(fbands{c})
end
Normalize_Clim(gcf)
savefig('Fig2b.fig')
close

%% Figure 2c. Euclidean distance?

%% Figure 4a. NA plots
splitmethod = {'amp','cos'};
for q = 1:2
    figure
    for c = 1:settings.nfreqs
        subplot(2,4,c)
        t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
        %         x2 = [t, fliplr(t)];
        %         inBetween = [mean(mean(allmeas{c}.naddttv.(splitmethod{q}).real(:,:,1,:),4),1), ...
        %             fliplr(mean(mean(allmeas{c}.naddttv.(splitmethod{q}).pseudo(:,:,1,:),4),1))];
        %         fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.3)
        hold on
        plot(t,mean(mean(allmeas{c}.naddttv.(splitmethod{q}).real(:,:,1,:),4),1),'b','LineWidth',1.5)
        plot(t,mean(mean(allmeas{c}.naddttv.(splitmethod{q}).pseudo(:,:,1,:),4),1),'b--','LineWidth',1.5)
        
        %         inBetween = [mean(mean(allmeas{c}.naddttv.(splitmethod{q}).real(:,:,2,:),4),1),...
        %             fliplr(mean(mean(allmeas{c}.naddttv.(splitmethod{q}).pseudo(:,:,2,:),4),1))];
        %         fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.3)
        plot(t,mean(mean(allmeas{c}.naddttv.(splitmethod{q}).pseudo(:,:,2,:),4),1),'r--','LineWidth',1.5)
        plot(t,mean(mean(allmeas{c}.naddttv.(splitmethod{q}).real(:,:,2,:),4),1),'r','LineWidth',1.5)
        title(fbands{c})
        legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'})
        xlabel('Time (s)')
        ylabelunits(settings)
    end
    savefig(['Fig4a_1_' num2str(q)])
    close
end

splitmethod = {'amp','cos'};
for q = 1:2
    figure
    for c = 1:settings.nfreqs
        subplot(2,4,c)
        t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
        %         x2 = [t, fliplr(t)];
        %         inBetween = [mean(mean(allmeas{c}.naddersp.(splitmethod{q}).real(:,:,1,:),4),1), ...
        %             fliplr(mean(mean(allmeas{c}.naddersp.(splitmethod{q}).pseudo(:,:,1,:),4),1))];
        %         fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.3)
        hold on
        plot(t,mean(mean(allmeas{c}.naddersp.(splitmethod{q}).real(:,:,1,:),4),1),'b','LineWidth',1.5)
        plot(t,mean(mean(allmeas{c}.naddersp.(splitmethod{q}).pseudo(:,:,1,:),4),1),'b--','LineWidth',1.5)
        
        %         inBetween = [mean(mean(allmeas{c}.naddersp.(splitmethod{q}).real(:,:,2,:),4),1),...
        %             fliplr(mean(mean(allmeas{c}.naddersp.(splitmethod{q}).pseudo(:,:,2,:),4),1))];
        %         fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.3)
        plot(t,mean(mean(allmeas{c}.naddersp.(splitmethod{q}).pseudo(:,:,2,:),4),1),'r--','LineWidth',1.5)
        plot(t,mean(mean(allmeas{c}.naddersp.(splitmethod{q}).real(:,:,2,:),4),1),'r','LineWidth',1.5)
        title(fbands{c})
        legend({'Real prestim low','Pseudo prestim low','Real prestim high','Pseudo prestim high'})
        xlabel('Time (s)')
        ylabelunits(settings)
        
    end
    savefig(['Fig4a_2_' num2str(q)])
    close
end

%% Figure 4b. Significance of NAindex for ERSP and TTV

for q = 1:2
    figure
    for c = 1:settings.nfreqs
        ax(c) = subplot(4,4,c);
        if strcmpi(settings.datatype,'MEG')
            ft_topoplot_vec(settings.layout,median(allmeas{c}.nattvindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ttv.stats{c}.mask),settings.datasetinfo.label);
        else
            topoplot(median(allmeas{c}.nattvindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ttv.stats{c}.mask),settings.layout);
        end
        title(fbands{c})
        
        ax(c+settings.nfreqs) = subplot(4,4,c+8);
        if strcmpi(settings.datatype,'MEG')
            ft_topoplot_vec(settings.layout,median(allmeas{c}.naerspindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ersp.stats{c}.mask),settings.datasetinfo.label);
        else
            topoplot(median(allmeas{c}.naerspindex.(splitmethod{q}),2).*...
                (alloutputs.([splitmethod{q} '_na']).ersp.stats{c}.mask),settings.layout);
        end
        title(fbands{c})
    end
    Normalize_Clim(gcf)
    colorbar
    savefig(['Fig4b_' num2str(q)])
    close
end



%% Figure 5: Resting state correlations



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

