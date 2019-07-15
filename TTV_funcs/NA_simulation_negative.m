% Simulation test for nonadditivity

mvals = 2.^[-6:6];
nvals = 2.^[-6:6];

for m = mvals
    for n = nvals
        datacalc = cell(1,48);
        
        parfor c = 1:48
            cfg = []; cfg.method = 'mix_oscifrac'; cfg.fsample = 500; cfg.trllen = 4;
            cfg.osci = struct; cfg.frac.ple = rand+0.5; cfg.frac.ampl = 1; cfg.frac.bpfreq = [0.5 50];
            %cfg.osci.s1.ampl = [ones(1,1000) -0.3*sin((1:1000)*pi/(1000))+1];
            cfg.numtrl = 128;
            cfg.noise.ampl = m + n*rand(1,128);
            for cc = 1:128
                spont = createFN(1.75/2,2000);
                spont = ft_preproc_lowpassfilter(spont,500,10,4);
                spont = NormOntoRange(spont,[1 2]);
                cfg.osci.s1.ampl{cc} = horz(spont)-[zeros(1,1000) rand*sin((1:1000)*pi/1000)];
            end
            sim{c} = ft_freqsimulation_swt(cfg);
            
            tmpcfg = []; tmpcfg.bpfilter = 'yes'; tmpcfg.bpfreq = [8 13]; tmpcfg.hilbert = 'complex';
            sim_bp{c} = ft_preprocessing(tmpcfg,sim{c});
            
            settings = struct;
            settings.units = 'prcchange';
            datacalc_tmp = Calc_sub(settings,sim_bp{c});
            datacalc{c} = datacalc_tmp{1};
        end
        all_datacalc{find(mvals == m,1),find(nvals == n,1),:} = datacalc;
    end
end

for c = 1:length(mvals)
    for cc = 1:length(nvals)
        all_datacalc{c,cc} = mergestructs(all_datacalc{c,cc});
    end
end

settings = parload('settings_camcan_1Hz.mat','settings');

datasetinfo = settings.datasetinfo;
datasetinfo.label = datasetinfo.label(1);

allstats_pt = cell(length(mvals),length(nvals));
allstats_ttv = allstats_pt;
for c = 1:length(mvals)
    opts.minnbchan = 0; opts.nrand = 1000;
    stats = cell(1,length(nvals));
    parfor cc = 1:length(nvals)
        stats_pt{cc} = EasyClusterCorrect({permute(squeeze(all_datacalc{c,cc}.naddersp.diff(1,:,1,:)),[3 2 1]),...
            permute(squeeze(all_datacalc{c,cc}.naddersp.diff(1,:,2,:)),[3 2 1])},datasetinfo,'ft_statfun_fast_signrank',opts);
        stats_ttv{cc} = EasyClusterCorrect({permute(all_datacalc{c,cc}.ttversp.real(1,:,:),[1 3 2]) 0.*permute(all_datacalc{c,cc}.ttversp.real(1,:,:),[1 3 2])},...
            datasetinfo,'ft_statfun_fast_signrank',opts)
    end
    allstats_pt(c,:) = stats_pt;
    allstats_ttv(c,:) = stats_ttv;
end

save('simulation_alloutputs.mat','allstats_pt','allstats_ttv','all_datacalc','-v7.3')

p = panel('no-manage-font');

pos = get(gcf,'position');

set(gcf,'position',[pos(1) pos(2) pos(3)*3 pos(4)*3]);

p.pack(13,13)

p.de.margin = [7 7 7 7];
p.marginleft = 18;

t = linspace(0,800,400);
warning('Hard-coded')
for c = 1:length(mvals)
    for cc = 1:length(nvals)
        p(c,cc).select()
        ylim1 = stdshade(t,squeeze(all_datacalc{c,cc}.naddersp.diff(1,:,1,:)),[0 0 1],0.15,1,'std');
        hold on
        ylim2 = stdshade(t,squeeze(all_datacalc{c,cc}.naddersp.diff(1,:,2,:)),[1 0 0],0.15,1,'std');
        set(gca,'XTickLabel',{},'XLim',[0 800],'YLim',[min(ylim1(1),ylim2(1)) max(ylim1(2),ylim2(2))])
                ytick = get(gca,'YTickLabel');
        if length(ytick) > 4
           set(gca,'YTickLabel',ytick([1:2:length(ytick)])) 
        end
        Plot_sigmask(gca,allstats_pt{c,cc}.prob < 0.05,'bar')
        set(gca,'XLim',[0 800])
        if c == 1
            title(['Noise TTV ' num2str(round(0.3*nvals(cc),3,'significant'))],'FontSize',10)
        end
        if cc == 1
            AddFigureLabel(gca,['SNR ' newline num2str(round(1/mvals(c),3,'significant'))],'middle_left','FontSize',10)
        
        end
        %         if c == length(mvals)
        %            xlabel('Time (ms)')
        %         end
        %         if cc == 1
        %            ylabel('% change from prestim')
        %         end
    end
end

set(gcf,'Color','w')

savefig('Simulation_negative_fig_pseudotrial.fig')
export_fig('Simulation_negative_fig_pseudotrial.png','-m4')



p = panel('no-manage-font');

pos = get(gcf,'position');

set(gcf,'position',[pos(1) pos(2) pos(3)*3 pos(4)*3]);

p.pack(13,13)

p.de.margin = [7 7 7 7];

p.marginleft = 18;

t = linspace(0,800,400);
warning('Hard-coded')
for c = 1:length(mvals)
    for cc = 1:length(nvals)
        p(c,cc).select()
        ylim1 = stdshade(t,squeeze(all_datacalc{c,cc}.ttversp.real(1,:,:)),[0 0 1],0.15,1,'std');
        set(gca,'XTickLabel',{},'XLim',[0 800],'YLim',ylim1)
        ytick = get(gca,'YTickLabel');
        if length(ytick) > 4
           set(gca,'YTickLabel',ytick([1:2:length(ytick)])) 
        end
        Plot_sigmask(gca,allstats_pt{c,cc}.prob < 0.05,'bar')
        set(gca,'XLim',[0 800])
        if c == 1
            title(['Noise TTV ' num2str(round(0.3*nvals(cc),3,'significant'))],'FontSize',10)
        end
        if cc == 1
            AddFigureLabel(gca,['SNR ' newline num2str(round(1/mvals(c),3,'significant'))],'middle_left','FontSize',10)
        end
        %         if c == length(mvals)
        %            xlabel('Time (ms)')
        %         end
        %         if cc == 1
        %            ylabel('% change from prestim')
        %         end
    end
end

set(gcf,'Color','w')

savefig('Simulation_negative_fig_ttv.fig')
export_fig('Simulation_negative_fig_ttv.png','-m4')


function [datacalc] = Calc_sub(settings,sim)

timefreq_data{1} = sim;

numbands = 1;

aucindex = 1:400;

datacalc = cell(1,1);
datacalc{1} = struct;

prestim_pseudo = 351:400; poststim_pseudo = 401:800; prestim_real = 951:1000; poststim_real = 1001:1400;

for q = 1:numbands
    nbchan = length(timefreq_data{q}.label);
    timefreq_data{q}.matrix = cat(3,timefreq_data{q}.trial{:});
    
    %% ERSP NA
    datacat = timefreq_data{q}.matrix;
    
    split_real = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_real,:)),2));
    split_pseudo = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_pseudo,:)),2));
    
    for c = 1:nbchan
        splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
        
        datacalc{q}.naddersp.raw.pseudo(c,:,1) = mean(abs(datacat(c,:,find(~splitindex))),3);
        datacalc{q}.naddersp.raw.pseudo(c,:,2) = mean(abs(datacat(c,:,find(splitindex))),3);
        
        datacalc{q}.naddersp.pseudo(c,:,1) = (mean(abs(datacat(c,poststim_pseudo,find(~splitindex))),3)...
            -mean(mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),3),2));
        datacalc{q}.naddersp.pseudo(c,:,2) = (mean(abs(datacat(c,poststim_pseudo,find(splitindex))),3)...
            -mean(mean(abs(datacat(c,prestim_pseudo,find(splitindex))),3),2));
        
        tmp = abs(datacat(c,poststim_pseudo,find(~splitindex)))...
            -mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),2);
        tmppseudo{1} = squeeze(trapz(tmp,2));
        
        tmp = abs(datacat(c,poststim_pseudo,find(splitindex)))...
            -mean(abs(datacat(c,prestim_pseudo,find(splitindex))),2);
        tmppseudo{2} = squeeze(trapz(tmp,2));
        
        switch settings.units
            case 'prcchange'
                datacalc{q}.naddersp.pseudo(c,:,1) = 100*datacalc{q}.naddersp.pseudo(c,:,1)./mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
                datacalc{q}.naddersp.pseudo(c,:,2) = 100*datacalc{q}.naddersp.pseudo(c,:,2)./mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
            case 'zscore'
                datacalc{q}.naddersp.pseudo(c,:,1) = zscore(datacalc{q}.naddersp.pseudo(c,:,1),0,2);
                datacalc{q}.naddersp.pseudo(c,:,2) = zscore(datacalc{q}.naddersp.pseudo(c,:,2),0,2);
            case 'log'
                datacalc{q}.naddersp.pseudo(c,:,1) = 10*log10(datacalc{q}.naddersp.pseudo(c,:,1));
                datacalc{q}.naddersp.pseudo(c,:,2) = 10*log10(datacalc{q}.naddersp.pseudo(c,:,2))
        end
        
        splitindex = split_real(c,:) > median(split_real(c,:));
        
        datacalc{q}.naddersp.raw.real(c,:,1) = mean(abs(datacat(c,:,find(~splitindex))),3);
        datacalc{q}.naddersp.raw.real(c,:,2) = mean(abs(datacat(c,:,find(splitindex))),3);
        
        datacalc{q}.naddersp.real(c,:,1) = (mean(abs(datacat(c,poststim_real,find(~splitindex))),3)...
            -mean(mean(abs(datacat(c,prestim_real,find(~splitindex))),3),2));
        datacalc{q}.naddersp.real(c,:,2) = (mean(abs(datacat(c,poststim_real,find(splitindex))),3)...
            -mean(mean(abs(datacat(c,prestim_real,find(splitindex))),3),2));
        
        datacalc{q}.naddersp.diff = datacalc{q}.naddersp.real- datacalc{q}.naddersp.pseudo;
        
        switch settings.units
            case 'prcchange'
                datacalc{q}.naddersp.real(c,:,1) = 100*datacalc{q}.naddersp.real(c,:,1)./mean(mean(abs(datacat(c,prestim_real,:)),3),2);
                datacalc{q}.naddersp.real(c,:,2) = 100*datacalc{q}.naddersp.real(c,:,2)./mean(mean(abs(datacat(c,prestim_real,:)),3),2);
            case 'zscore'
                datacalc{q}.naddersp.real(c,:,1) = zscore(datacalc{q}.naddersp.real(c,:,1),0,2);
                datacalc{q}.naddersp.real(c,:,2) = zscore(datacalc{q}.naddersp.real(c,:,2),0,2);
            case 'log'
                datacalc{q}.naddersp.real(c,:,1) = 10*log10(datacalc{q}.naddersp.real(c,:,1));
                datacalc{q}.naddersp.real(c,:,2) = 10*log10(datacalc{q}.naddersp.real(c,:,2));
        end
    end
    %% TTV of ERSP
    
    datacalc{q}.raw.ttversp(:,:) = std(abs(datacat),[],3);
    datacalc{q}.ttversp.pseudo(:,:) = (std(abs(datacat(:,poststim_pseudo,:)),[],3)...
        -mean(std(abs(datacat(:,prestim_pseudo,:)),[],3),2));
    datacalc{q}.ttversp.real(:,:) = (std(abs(datacat(:,poststim_real,:)),[],3)...
        -mean(std(abs(datacat(:,prestim_real,:)),[],3),2));
    
    switch settings.units
        case 'prcchange'
            datacalc{q}.ttversp.real(:,:) = 100*datacalc{q}.ttversp.real(:,:)./mean(datacalc{q}.raw.ttversp(:,prestim_real),2);
            datacalc{q}.ttversp.pseudo(:,:) = 100*datacalc{q}.ttversp.pseudo(:,:)./mean(datacalc{q}.raw.ttversp(:,prestim_pseudo),2);
        case 'zscore'
            datacalc{q}.ttversp.real(:,:) = zscore(datacalc{q}.ttversp.real(:,:),0,2);
            datacalc{q}.ttversp.pseudo(:,:) = zscore(datacalc{q}.ttversp.pseudo(:,:),0,2);
        case 'log'
            datacalc{q}.ttversp.real(:,:) = 10*log10(datacalc{q}.ttversp.real(:,:));
            datacalc{q}.ttversp.pseudo(:,:) = 10*log10(datacalc{q}.ttversp.pseudo(:,:));
    end
    
    
end
end
