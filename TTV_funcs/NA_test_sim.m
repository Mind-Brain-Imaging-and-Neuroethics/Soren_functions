% Simulation test for nonadditivity


parfor i = 1:50
    cfg = []; cfg.method = 'mix_oscifrac'; cfg.fsample = 500; cfg.trllen = 4;
    cfg.osci = struct; cfg.frac.ple = rand+0.3; cfg.frac.ampl = 1; cfg.frac.bpfreq = [0.5 50];
    simtrl = cell(1,128);
    for c = 1:128
        %pseudoampl = ft_preproc_bandpassfilter(randn(1,1000),500,[0.1 1]);
        %realampl = linspace(1,0.5,1000)*ft_preproc_bandpassfilter(randn(1,1000),500,[0.1 1]);
        cfg.osci.s1.ampl = [ones(1,1000) -0.3*sin((1:1000)*pi/(2000))+1];
        cfg.noise.ampl = 5*rand+1;
        cfg.numtrl = 1;
        simtrl{c} = ft_freqsimulation(cfg);
    end
    tmpcfg = [];
    sim{i} = ft_appenddata(tmpcfg,simtrl{:});
    
    cfg = []; cfg.bpfilter = 'yes'; cfg.bpfreq = [8 13]; cfg.hilbert = 'complex';
    sim_bp{i} = ft_preprocessing(cfg,sim{i});
    
    settings = struct;
    settings.units = 'prcchange';
    datacalc{i} = Calc_sub(settings,sim_bp{i});
    datacalc{i} = datacalc{i}{1};
end

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
    
    SD_all = zeros(length(timefreq_data{q}.label),size(timefreq_data{q}.matrix,2));
    for c = 1:length(timefreq_data{q}.label)
        SD_all(c,:) = std(real(timefreq_data{q}.matrix(c,:,:)),[],3);
    end
    datacalc{q}.raw.sd(:,:) = SD_all;
    
    % normalize by 50ms prestim ifor both real and pseudotrial
    datacalc{q}.ttv.pseudo(:,:) = (SD_all(:,poststim_pseudo)-...
        mean(SD_all(:,prestim_pseudo),2));
    datacalc{q}.ttv.real(:,:) = (SD_all(:,poststim_real)-...
        mean(SD_all(:,prestim_real),2));
    switch settings.units
        case 'prcchange'
            datacalc{q}.ttv.pseudo(:,:) = 100*datacalc{q}.ttv.pseudo(:,:)./...
                mean(SD_all(:,prestim_pseudo),2);
            datacalc{q}.ttv.real(:,:) = 100*datacalc{q}.ttv.real(:,:)./...
                mean(SD_all(:,prestim_real),2);
        case 'zscore'
            datacalc{q}.ttv.pseudo(:,:) = zscore(datacalc{q}.ttv.pseudo(:,:),0,2);
            datacalc{q}.ttv.real(:,:) = zscore(datacalc{q}.ttv.real(:,:),0,2);
        case 'log'
            datacalc{q}.ttv.pseudo(:,:) = 10*log10(datacalc{q}.ttv.pseudo(:,:));
            datacalc{q}.ttv.real(:,:) = 10*log10(datacalc{q}.ttv.real(:,:));
    end
    
    
    %% ERSP and ITC
    datacat = timefreq_data{q}.matrix;
    
    datacalc{q}.raw.ersp(:,:) = mean(abs(timefreq_data{q}.matrix),3);
    datacalc{q}.raw.itc(:,:) = abs(mean(datacat./abs(datacat),3));
    
    datacalc{q}.ersp.pseudo(:,:) = (mean(abs(datacat(:,poststim_pseudo,:)),3)...
        -mean(mean(abs(datacat(:,prestim_pseudo,:)),3),2));
    datacalc{q}.ersp.real(:,:) = (mean(abs(datacat(:,poststim_real,:)),3)...
        -mean(mean(abs(datacat(:,prestim_real,:)),3),2));
    
    switch settings.units
        case 'prcchange'
            datacalc{q}.ersp.pseudo(:,:) = 100*datacalc{q}.ersp.pseudo(:,:)./mean(mean(abs(datacat(:,prestim_pseudo,:)),3),2);
            datacalc{q}.ersp.real(:,:) = 100*datacalc{q}.ersp.real(:,:)./mean(mean(abs(datacat(:,prestim_real,:)),3),2);
        case 'zscore'
            datacalc{q}.ersp.pseudo(:,:) = zscore(datacalc{q}.ersp.pseudo(:,:),0,2);
            datacalc{q}.ersp.real(:,:) = zscore(datacalc{q}.ersp.real(:,:),0,2);
        case 'log'
            datacalc{q}.ersp.pseudo(:,:) = 10*log10(datacalc{q}.ersp.pseudo(:,:));
            datacalc{q}.ersp.real(:,:) = 10*log10(datacalc{q}.ersp.real(:,:));
    end
    
    datacalc{q}.itc.real(:,:) = (datacalc{q}.raw.itc(:,poststim_real)-mean(datacalc{q}.raw.itc(:,prestim_real),2));
    datacalc{q}.itc.pseudo(:,:) = (datacalc{q}.raw.itc(:,poststim_pseudo)-mean(datacalc{q}.raw.itc(:,prestim_pseudo),2))./...
        mean(datacalc{q}.raw.itc(:,prestim_pseudo),2);
    
    switch settings.units
        case 'prcchange'
            datacalc{q}.itc.real(:,:) = 100*datacalc{q}.itc.real(:,:)./mean(datacalc{q}.raw.itc(:,prestim_real),2);
            datacalc{q}.itc.pseudo(:,:) = 100*datacalc{q}.itc.pseudo(:,:)./mean(datacalc{q}.raw.itc(:,prestim_pseudo),2);
        case 'zscore'
            datacalc{q}.itc.real(:,:) = zscore(datacalc{q}.itc.real(:,:),0,2);
            datacalc{q}.itc.pseudo(:,:) = zscore(datacalc{q}.itc.pseudo(:,:),0,2);
        case 'log'
            datacalc{q}.itc.real(:,:) = 10*log10(datacalc{q}.itc.real(:,:));
            datacalc{q}.itc.pseudo(:,:) = 10*log10(datacalc{q}.itc.pseudo(:,:));
    end
    
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
        
        %         tmp = abs(datacat(c,poststim_pseudo,find(~splitindex)))...
        %             -mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),2);
        %         tmpreal{1}(:) = squeeze(trapz(tmp,2));
        %
        %         tmp = abs(datacat(c,poststim_pseudo,find(splitindex)))...
        %             -mean(abs(datacat(c,prestim_pseudo,find(splitindex))),2);
        %         tmpreal{2}(:) = squeeze(trapz(tmp,2));
        
        %        [~,~,~,stats] = ttest2(tmpreal{2}-tmppseudo{2},tmpreal{1}-tmppseudo{1});
        %        datacalc{q}.t(c) = stats.t;
        
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
    
    %% ERP nonadditivity
    datacalc{q}.erp(:,:) = mean(real(timefreq_data{q}.matrix),3);
    
    split_real = squeeze(mean(real(timefreq_data{q}.matrix(:,prestim_real,:)),2));
    split_pseudo = squeeze(mean(real(timefreq_data{q}.matrix(:,prestim_pseudo,:)),2));
    
    %datacat = timefreq_data{q}.matrix;
    
    for c = 1:nbchan
        splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
        
        datacalc{q}.nadderp.raw.pseudo(c,:,1) = mean(real(datacat(c,:,find(~splitindex))),3);
        datacalc{q}.nadderp.raw.pseudo(c,:,2) = mean(real(datacat(c,:,find(splitindex))),3);
        
        
        datacalc{q}.nadderp.pseudo(c,:,1) = (mean(real(datacat(c,poststim_pseudo,find(~splitindex))),3)...
            -mean(mean(real(datacat(c,prestim_pseudo,find(~splitindex))),3),2));
        datacalc{q}.nadderp.pseudo(c,:,2) = (mean(real(datacat(c,poststim_pseudo,find(splitindex))),3)...
            -mean(mean(real(datacat(c,prestim_pseudo,find(splitindex))),3),2));
        
        %             switch settings.units
        %                 case 'prcchange'
        %                     datacalc{q}.nadderp.pseudo(c,:,1) = 100*datacalc{q}.nadderp.pseudo(c,:,1)./mean(mean(real(datacat(c,prestim_pseudo,:)),3),2);
        %                     datacalc{q}.nadderp.pseudo(c,:,2) = 100*datacalc{q}.nadderp.pseudo(c,:,2)./mean(mean(real(datacat(c,prestim_pseudo,:)),3),2);
        %                 case 'zscore'
        %                     datacalc{q}.nadderp.pseudo(c,:,1) = zscore(datacalc{q}.nadderp.pseudo(c,:,1),0,2);
        %                     datacalc{q}.nadderp.pseudo(c,:,2) = zscore(datacalc{q}.nadderp.pseudo(c,:,2),0,2);
        %                 case 'log'
        %                     datacalc{q}.nadderp.pseudo(c,:,1) = 10*log10(datacalc{q}.nadderp.pseudo(c,:,1));
        %                     datacalc{q}.nadderp.pseudo(c,:,2) = 10*log10(datacalc{q}.nadderp.pseudo(c,:,2));
        %             end
        
        splitindex = split_real(c,:) > median(split_real(c,:));
        
        datacalc{q}.nadderp.raw.real(c,:,1) = mean(real(datacat(c,:,find(~splitindex))),3); %NOT ACTUALLY THE RAW ERSP - used for plotting later
        datacalc{q}.nadderp.raw.real(c,:,2) = mean(real(datacat(c,:,find(splitindex))),3);
        
        
        datacalc{q}.nadderp.real(c,:,1) = (mean(real(datacat(c,poststim_real,find(~splitindex))),3)...
            -mean(mean(real(datacat(c,prestim_real,find(~splitindex))),3),2));
        datacalc{q}.nadderp.real(c,:,2) = (mean(real(datacat(c,poststim_real,find(splitindex))),3)...
            -mean(mean(real(datacat(c,prestim_real,find(splitindex))),3),2));
        
        %             switch settings.units
        %                 case 'prcchange'
        %                     datacalc{q}.nadderp.real(c,:,1) = 100*datacalc{q}.nadderp.real(c,:,1)./mean(mean(real(datacat(c,prestim_real,:)),3),2);
        %                     datacalc{q}.nadderp.real(c,:,2) = 100*datacalc{q}.nadderp.real(c,:,2)./mean(mean(real(datacat(c,prestim_real,:)),3),2);
        %                 case 'zscore'
        %                     datacalc{q}.nadderp.real(c,:,1) = zscore(datacalc{q}.nadderp.real(c,:,1),0,2);
        %                     datacalc{q}.nadderp.real(c,:,2) = zscore(datacalc{q}.nadderp.real(c,:,2),0,2);
        %                 case 'log'
        %                     datacalc{q}.nadderp.real(c,:,1) = 10*log10(datacalc{q}.nadderp.real(c,:,1));
        %                     datacalc{q}.nadderp.real(c,:,2) = 10*log10(datacalc{q}.nadderp.real(c,:,2))
        %             end
        
        datacalc{q}.nadderp.diff = datacalc{q}.nadderp.real-datacalc{q}.nadderp.pseudo;
    end
end
end
