function camcan_oscifrac_figure2

addpath(genpath('/project/def-gnorthof/sorenwt/MATLAB'))
addpath('/project/def-gnorthof/sorenwt/fieldtrip-master')
ft_defaults
addpath([toolboxdir('signal'),'/signal'])

load('settings_camcan_1Hz.mat')

cd /scratch/sorenwt/camcan/Preprocessed/Task/Epoched/

files = dir('*tf.mat');

m = matfile(files(1).name);

fields = {'mixd','osci','frac'};

for c = 1:length(fields)
    meandata.(fields{c}) = m.(fields{c});
end

%% Reading in data
donefiles = [];
for i = 1:length(files)
    try
        %if ~ismember({files(i).name},donefiles)
            m = matfile(files(i).name);
            for c = 1:3
                tmp = m.(fields{c});
                meandata.(fields{c}).fourierspctrm(i,:,:,:) = nanmean(tmp.fourierspctrm,1);
            end
            donefiles = [donefiles {files(i).name}];
        %end
    catch
    end
end

for c = 1:length(fields)
   meandata.(fields{c}).fourierspctrm = permute(meandata.(fields{c}).fourierspctrm,[1 3 2 4]); 
   meandata.(fields{c}).fourierspctrm_dimord = meandata.(fields{c}).dimord;
   meandata.(fields{c}).grad = settings.datasetinfo.grad;
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

post = cell(1,size(meanpost.frac.fourierspctrm,1));
bl = post;

parfor i = 1:size(meanpost.frac.fourierspctrm,1)
    for ii = 1:size(meanpost.frac.fourierspctrm,2)
        for iii = 1:size(meanpost.frac.fourierspctrm,4)
            tmp = log10(meanpost.frac.freq);
            pow = log10(squeeze(meanpost.frac.fourierspctrm(i,ii,:,iii)));
            lintmp = vert(linspace(tmp(1),tmp(end),length(tmp)));
            pow = interp1(tmp,pow,lintmp);
            p = polyfit(lintmp,pow,1);
            post{i}(1,ii,iii) = -p(1);
            
            tmp = log10(meanbl.frac.freq);
            pow = log10(squeeze(meanbl.frac.fourierspctrm(i,ii,:,iii)));
            lintmp = vert(linspace(tmp(1),tmp(end),length(tmp)));
            pow = interp1(tmp,pow,lintmp);
            p = polyfit(lintmp,pow,1);
            bl{i}(1,ii,iii) = -p(1);
        end
    end
end

pledata.post = cat(1,post{:});
pledata.bl = cat(1,bl{:});

for c = 1:3
    cfg = []; cfg.channel = {'MEG'}; cfg.avgoverchan = 'yes'; cfg.latency = [0 0.75];
    cfg.frequency = 'all'; cfg.method = 'montecarlo'; cfg.statistic = 'ft_statfun_actvsblT';
    cfg.correctm = 'cluster'; cfg.clusteralpha = 0.05; cfg.clusterstatistic = 'maxsum';
    cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025; cfg.numrandomization = 2000;
    cfg.parameter = 'fourierspctrm';
    
    ntrials = size(meanpost.(fields{c}).fourierspctrm,1);
    design  = zeros(2,2*ntrials);
    design(1,1:ntrials) = 1;
    design(1,ntrials+1:2*ntrials) = 2;
    design(2,1:ntrials) = [1:ntrials];
    design(2,ntrials+1:2*ntrials) = [1:ntrials];
    
    cfg.design = design;
    cfg.ivar = 1;
    cfg.uvar = 2;
    
    cfg.parpool = 48;
    
    meanbl.(fields{c}).time = meanpost.(fields{c}).time;
    
    stats{c} = ft_freqstatistics(cfg,meanpost.(fields{c}),meanbl.(fields{c}));
end

% Cluster stats for PLE and broadband

datasetinfo = settings.datasetinfo;

opts = []; opts.nrand = 2000; opts.parpool = 48; opts.minnbchan = 0;
opts.external = cfg.design;

stats_ple = EasyClusterCorrect({permute(pledata.post,[2 1 3]) repmat(mean(pledata.bl,3)',1,1,38)},...
    datasetinfo,'ft_statfun_signrank',opts);

stats_bb = EasyClusterCorrect({permute(squeeze(mean(meanpost.frac.fourierspctrm,3)),[2 1 3]) repmat(mean(squeeze(mean(meanbl.frac.fourierspctrm,3)),3)',1,1,38)},...
    datasetinfo,'ft_statfun_signrank',opts);



for c = find(extractfield(stats{1}.posclusters,'prob') < 0.05)
    statmask = squeeze(stats{1}.posclusterslabelmat==c);
    % for ple and broadband fractal, include those time points where a
    % number of frequencies greater than the median value in the cluster
    % are significant
    sum_statmask = sum(statmask,1);
    %statmask_time = sum_statmask > median(sum_statmask(find(sum_statmask>0))); 
    statmask_time = sum_statmask > 0;
    meanmix(:,c) = sum(sum(permute(permute(squeeze(mean(meanpost.mixd.fourierspctrm,2)...
        -mean(mean(meanbl.mixd.fourierspctrm,4),2)),[2 3 1]).*statmask,[3 1 2]),2),3); % baseline correct, sum over time and freq for significant cluster
    meanosci(:,c) = sum(sum(permute(permute(squeeze(mean(meanpost.osci.fourierspctrm,2)...
        -mean(mean(meanbl.osci.fourierspctrm,4),2)),[2 3 1]).*statmask,[3 1 2]),2),3);
    meanfrac(:,c) = mean(squeeze(mean(mean(meanpost.frac.fourierspctrm,2),3)...
        -mean(mean(mean(meanbl.frac.fourierspctrm,4),3),2)).*statmask_time,2); %use only the time statmask here, not freq - broadband power
    meanple(:,c) = sum(squeeze(mean(pledata.post,2)-mean(mean(pledata.bl,2),3)).*statmask_time,2);
end

for c = find(extractfield(stats{1}.negclusters,'prob') < 0.05)
    statmask = squeeze(stats{1}.negclusterslabelmat==c);
    % for ple and broadband fractal, include those time points where a
    % number of frequencies greater than the median value in the cluster
    % are significant
    sum_statmask = sum(statmask,1);
    %statmask_time = sum_statmask > median(sum_statmask(find(sum_statmask>0))); 
    statmask_time = sum_statmask > 0;
    
    newindx = c+length(find(extractfield(stats{1}.posclusters,'prob') < 0.05));
    % for ple and broadband fractal, include those time points where a
    % number of frequencies greater than the median value in the cluster
    % are significant
    sum_statmask = sum(statmask,1);
    statmask_time = sum_statmask > median(sum_statmask(find(sum_statmask>0))); 
    meanmix(:,newindx) = sum(sum(permute(permute(squeeze(mean(meanpost.mixd.fourierspctrm,2)...
        -mean(mean(meanbl.mixd.fourierspctrm,4),2)),[2 3 1]).*statmask,[3 1 2]),2),3); % baseline correct, sum over time and freq for significant cluster
    meanosci(:,newindx) = sum(sum(permute(permute(squeeze(mean(meanpost.osci.fourierspctrm,2)...
        -mean(mean(meanbl.osci.fourierspctrm,4),2)),[2 3 1]).*statmask,[3 1 2]),2),3);
    meanfrac(:,newindx) = mean(squeeze(mean(mean(meanpost.frac.fourierspctrm,2),3)...
        -mean(mean(mean(meanbl.frac.fourierspctrm,4),3),2)).*statmask_time,2); %use only the time statmask here, not freq - broadband power
    meanple(:,newindx) = sum(squeeze(mean(pledata.post,2)-mean(mean(pledata.bl,2),3)).*statmask_time,2);
end

%% Regression models

for i = 1:size(meanmix,2)
   tbl{i} = array2table([meanmix(:,i) meanosci(:,i) meanple(:,i) meanfrac(:,i)],'VariableNames',{'Mixed_power','Osci_power','PLE','Frac_BB'});
   tbl{i}{:,:} = zscore(tbl{i}{:,:},[],1);
   mdl{i} = fitlm(tbl{i},'Mixed_power~Osci_power+PLE+Frac_BB');
   % partial correlation stuff
end

%% Plotting the figure

p = panel('no-manage-font');

p.pack('v',{50 50})

p(1).pack('h',{1/4 1/4 1/4 1/4})

for c = 1:3
    p(1,c).select()
    
    plotdata = meandata.(fields{c});
    plotdata.fourierspctrm(find(plotdata.fourierspctrm < 0)) = min(min(min(min(plotdata.fourierspctrm(find(plotdata.fourierspctrm > 0))))));
        plotdata.fourierspctrm(:,1,:,:) = mean(plotdata.fourierspctrm,2); % make the first channel the mean so you can use plotting mask
    plotdata.mask = logical(cat(3,zeros(1,size(stats{c}.mask,2),length(plotdata.time)-size(stats{c}.mask,3)),stats{c}.mask));
    plotdata.fourierspctrm_dimord = 'rpt_chan_freq_time';
    plotdata.dimord = 'rpt_chan_freq_time';
    
    cfg = []; cfg.parameter = 'fourierspctrm'; cfg.baseline = [-Inf 0]; cfg.baselinetype = 'db';
    cfg.maskparameter = 'mask'; cfg.maskstyle = 'outline';
    cfg.layout = 'neuromag306mag.lay'; cfg.channel = 1; cfg.latency = [0 0.75];
    cfg.interactive = 'no';
    ft_singleplotTFR(cfg,plotdata);
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    set(gca,'XLim',[-0.2 0.75])
    hold on
    line([0 0],get(gca,'YLim'),'Color',[0 0 0],'LineWidth',2)
    title(fields{c},'FontSize',14)
end

p(1,4).pack('v',{50 50})

p(1,4,1).select()
stdshade(meanpost.mixd.time,squeeze(mean(pledata.post,2)),'k',0.15,2,'sem')
Plot_sigmask(gca,stats_ple.mask,'cmapline')
FixAxes(gca,14)
xlabel('Time (s)')
ylabel('PLE')
set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)])

p(1,4,2).select()
stdshade(meanpost.mixd.time,squeeze(mean(mean(meanpost.frac.fourierspctrm,2),3)),'k',0.15,2,'sem')
Plot_sigmask(gca,stats_bb.mask,'cmapline')
FixAxes(gca,14)
xlabel('Time (s)')
ylabel('Fractal broadband power')
set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)])



p(2).pack('h',repmat({1/length(mdl)},1,length(mdl)));

for i = 1:length(mdl)
   p(2,i).select()
   scatter(1:3,mdl{i}.Coefficients.Estimate(2:end),36,[0 0 1],'x','LineWidth',2)
   hold on
   %scatter(1:3,partr(i,:),36,[1 0 0],'o','LineWidth',2)
   er = errorbar(1:3,mdl{i}.Coefficients.Estimate(2:end),mdl{i}.Coefficients.SE(2:end)*1.96,...
       'LineStyle','none','LineWidth',2,'Color',[0 0 1],'HandleVisibility','off');
   xl = get(gca,'XLim');
   line(xl+[-0.1 0.1],[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
   set(gca,'XLim',xl + [-0.1 0.1],'XTickLabel',{'Oscilatory Power','PLE','Fractal Broadband'})
   %legend({'Regression Coefficient','Partial R^2'})
   %ylabel('Regression Coefficient')
   FixAxes(gca,14)
   fix_xticklabels(gca,0.1,{'FontSize',14}) 
   ylabel('Coefficient')
end

p(1).de.marginleft = 30;
p.marginright = 20;
p.marginleft = 18;
p.margintop = 8;



