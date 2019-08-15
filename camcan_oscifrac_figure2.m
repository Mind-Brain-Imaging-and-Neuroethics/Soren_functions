function camcan_oscifrac_figure2

addpath(genpath('/project/def-gnorthof/sorenwt/MATLAB'))
addpath('/project/def-gnorthof/sorenwt/fieldtrip-master')
ft_defaults
addpath([toolboxdir('signal'),'/signal'])

cd /scratch/sorenwt/camcan/Preprocessed/Task/Epoched/

files = dir('*tf.mat');

m = matfile(files(1).name);

fields = {'mixd','osci','frac'};

for c = 1:length(fields)
    meandata.(fields{c}) = m.(fields{c});
end

%% Reading in data
for i = 1:length(files)
    m = matfile(files(i).name);
    for c = 1:3
        tmp = m.(fields{c});
        meandata.(fields{c}).fourierspctrm(i,:,:,:) = nanmean(tmp.fourierspctrm,1);
    end
end


%% Analysis and statistics

for c = 1:3
    post = find(meandata.(fields{c}).time >= 0);
    meanpost.(fields{c}) = meandata.(fields{c});
    meanpost.(fields{c}).fourierspctrm = meandata.(fields{c}).fourierspctrm(:,:,:,post);
    meanpost.(fields{c}).time = meandata.(fields{c}).time(post);
    
    bl = find(meandata.(fields{c}).time < 0);
    bl = bl((end-length(post)+1):end);
    meanbl.(fields{c}) = meandata.(fields{c});
    meanbl.(fields{c}).fourierspctrm = meandata.(fields{c}).fourierspctrm(:,:,:,bl);
    meanbl.(fields{c}).time = meandata.(fields{c}).time(bl);
end

for i = 1:size(meanpost.frac.fourierspctrm,1)
    for ii = 1:size(meanpost.frac.fourierspctrm,2)
        for iii = 1:size(meanpost.frac.fourierspctrm,4)
            tmp = log10(meanpost.frac.foi);
            pow = log10(squeeze(meanpost.frac.fourierspctrm(i,ii,:,iii)));
            lintmp = linspace(tmp(1),tmp(2),length(tmp));
            pow = interp1(tmp,pow,lintmp);
            p = polyfit(lintmp,pow,1);
            meanple.post(i,ii,iii) = -p(1);
            
            tmp = log10(meanbl.frac.foi);
            pow = log10(squeeze(meanbl.frac.fourierspctrm(i,ii,:,iii)));
            lintmp = linspace(tmp(1),tmp(2),length(tmp));
            pow = interp1(tmp,pow,lintmp);
            p = polyfit(lintmp,pow,1);
            meanple.bl(i,ii,iii) = -p(1);
        end
    end
end

for c = 1:3
    cfg = []; cfg.channel = {'MEG'}; cfg.avgoverchan = 'yes'; cfg.latency = [0 0.75];
    cfg.frequency = 'all'; cfg.method = 'montecarlo'; cfg.statistic = 'ft_statfun_actvsblT';
    cfg.correctm = 'cluster'; cfg.clusteralpha = 0.05; cfg.clusterstatistic = 'maxsum';
    cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025; cfg.numrandomization = 1000;
    
    ntrials = size(meanpost.(fields{c}).fourierspctrm,1);
    design  = zeros(2,2*ntrials);
    design(1,1:ntrials) = 1;
    design(1,ntrials+1:2*ntrials) = 2;
    design(2,1:ntrials) = [1:ntrials];
    design(2,ntrials+1:2*ntrials) = [1:ntrials];
    
    cfg.design = design;
    cfg.ivar = 1;
    cfg.uvar = 2;
    
    stats{c} = ft_freqstatistics(cfg,meanpost.(fields{c}),meanbl.(fields{c}));
end

for c = 1:length(unique(stats{1}.posclusterslabelmat))
    statmask = (stats{1}.posclusterslabelmat==c);
    % for ple and broadband fractal, include those time points where a
    % number of frequencies greater than the median value in the cluster
    % are significant
    sum_statmask = sum(statmask,1);
    statmask_time = sum_statmask > median(sum_statmask(find(sum_statmask>0))); 
    meanmix(:,c) = sum(sum((mean(meanpost.mixd.fourierspctrm,2)...
        -mean(mean(meanbl.mixd.fourierspctrm,4),2)).*statmask,3),4); % baseline correct, sum over time and freq for significant cluster
    meanosci(:,c) = sum(sum((mean(meanpost.osci.fourierspctrm,2)...
        -mean(mean(meanbl.osci.fourierspctrm,4),2)).*statmask,3),4);
    meanfrac(:,c) = mean((mean(mean(meanpost.frac.fourierspctrm,2),3)...
        -mean(mean(mean(meanbl.frac.fourierspctrm,4),3),2)).*statmask_time,4); %use only the time statmask here, not freq - broadband power
    meanple(:,c) = sum((mean(meanple.post,2)-mean(mean(meanple.bl,2),3)).*statmask_time,4);
end

for c = 1:length(unique(stats{1}.negclusterslabelmat))
    statmask = (stats{1}.negclusterslabelmat==c);
    % for ple and broadband fractal, include those time points where a
    % number of frequencies greater than the median value in the cluster
    % are significant
    sum_statmask = sum(statmask,1);
    statmask_time = sum_statmask > median(sum_statmask(find(sum_statmask>0))); 
    
    newindx = c+length(unique(stats{1}.posclusterslabelmat));
    meanmix(:,newindx) = sum(sum((mean(meanpost.mixd.fourierspctrm,2)...
        -mean(mean(meanbl.mixd.fourierspctrm,4),2)).*statmask,3),4); % baseline correct, sum over time and freq for significant cluster
    meanosci(:,newindx) = sum(sum((mean(meanpost.osci.fourierspctrm,2)...
        -mean(mean(meanbl.osci.fourierspctrm,4),2)).*statmask,3),4);
    meanfrac(:,newindx) = mean((mean(mean(meanpost.frac.fourierspctrm,2),3)...
        -mean(mean(mean(meanbl.frac.fourierspctrm,4),3),2)).*statmask_time,4); %use only the time statmask here, not freq - broadband power
    meanple(:,newindx) = sum((mean(meanple.post,2)-mean(mean(meanple.bl,2),3)).*statmask_time,4);
end

%% Regression models

for i = 1:size(meanmix,2)
   tbl{i} = array2table([meanmix(:,i) meanosci(:,i) meanple(:,i) meanfrac(:,i)],'VariableNames',{'Mixed_power','Osci_power','PLE','Frac_BB'});
   tbl{i}{:,:} = zscore(tbl{i}{:,:},[],1);
   mdl{i} = fitlm(tbl,'Mixed_power~Osci_power+PLE+Frac_BB');
   % partial correlation stuff
end

%% Plotting the figure

p = panel('no-manage-font');

p.pack('v',{50 50})

p(1).pack('h',{1/3 1/3 1/3})

for c = 1:3
    p(1,c).select()
    
    plotdata = meandata.(fields{c});
    plotdata.fourierspctrm(:,1,:,:) = mean(plotdata.fourierspctrm,2); % make the first channel the mean so you can use plotting mask
    plotdata.mask = stats{c}.mask;
    
    cfg = []; cfg.parameter = 'fourierspctrm'; cfg.baseline = [-Inf 0]; cfg.baselinetype = 'db';
    cfg.maskparameter = 'mask'; cfg.maskstyle = 'saturation';
    cfg.layout = 'neuromag306mag.lay'; cfg.channel = 1; cfg.latency = [0 0.75];
    cfg.interactive = 'no';
    ft_singleplotTFR(cfg,plotdata);
end

p(2).pack('h',repmat({1/length(mdl)},1,length(mdl)));

for c = 1:length(mdl)
   p(2,c).select()
   scatter(1:3,mdl{i}.Coefficients.Estimate(2:end),36,[0 0 1],'x','LineWidth',2)
   hold on
   %scatter(1:3,partr(i,:),36,[1 0 0],'o','LineWidth',2)
   er = errorbar(1:3,mdl{i}.Coefficients.Estimate(2:end),mdl{i}.Coefficients.SE(2:end)*1.96,...
       'LineStyle','none','LineWidth',2,'Color',[0 0 1],'HandleVisibility','off');
   xl = get(gca,'XLim');
   line(xl+[-0.1 0.1],[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
   set(gca,'XLim',xl + [-0.1 0.1],'XTickLabel',{'Osci_power','PLE','Frac BB'})
   %legend({'Regression Coefficient','Partial R^2'})
   %ylabel('Regression Coefficient')
   FixAxes(gca,14)
   fix_xticklabels(gca,0.1,{'FontSize',14}) 
end




