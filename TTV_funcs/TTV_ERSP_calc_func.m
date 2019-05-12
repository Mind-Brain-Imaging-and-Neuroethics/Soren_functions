function [settings] = TTV_ERSP_calc_func(settings)

numbands = settings.nfreqs;

aucindex = settings.aucindex;

if isempty(gcp('nocreate'))
    parpool(numbands)
end

files = dir(['*' settings.datasetname '_ft.mat_timefreq_filtered.mat']);

allmeas = cell(1,settings.nfreqs);
for c = 1:settings.nfreqs
    allmeas{c} = struct;
end

if strcmpi(settings.tfparams.method,'hilbert')
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

for i = 1:length(files)
    disp(['Processing subject ' num2str(i)])
    timefreq_data = parload(files(i).name,'timefreq_data');
    nbchan = size(timefreq_data{1}.trial{1},1);
    for q = 1:numbands
        
        timefreq_data{q}.matrix = cat(3,timefreq_data{q}.trial{:});
        
        SD_all = zeros(length(timefreq_data{q}.label),size(timefreq_data{q}.matrix,2));
        for c = 1:length(timefreq_data{q}.label)
            SD_all(c,:) = std(real(timefreq_data{q}.matrix(c,:,:)),[],3);
            %TTV_real(c,:) = std(read(timefreq_data{q}.matrix(c,191:380,:)),[],3);
        end
        allmeas{q}.raw.sd(:,:,i) = SD_all;
        
        % normalize by 50ms prestim ifor both real and pseudotrial
        allmeas{q}.ttv.pseudo(:,:,i) = (SD_all(:,poststim_pseudo)-...
            mean(SD_all(:,prestim_pseudo),2));
        allmeas{q}.ttv.real(:,:,i) = (SD_all(:,poststim_real)-...
            mean(SD_all(:,prestim_real),2));
        switch settings.units
            case 'prcchange'
                allmeas{q}.ttv.pseudo(:,:,i) = 100*allmeas{q}.ttv.pseudo(:,:,i)./...
                    mean(SD_all(:,prestim_pseudo),2);
                allmeas{q}.ttv.real(:,:,i) = 100*allmeas{q}.ttv.real(:,:,i)./...
                    mean(SD_all(:,prestim_real),2);
            case 'zscore'
                allmeas{q}.ttv.pseudo(:,:,i) = zscore(allmeas{q}.ttv.pseudo(:,:,i),0,2);
                allmeas{q}.ttv.real(:,:,i) = zscore(allmeas{q}.ttv.real(:,:,i),0,2);
            case 'log'
                allmeas{q}.ttv.pseudo(:,:,i) = 10*log10(allmeas{q}.ttv.pseudo(:,:,i));
                allmeas{q}.ttv.real(:,:,i) = 10*log10(allmeas{q}.ttv.real(:,:,i));
        end
        
        %% Nonadditivity TTV - amplitude based
        
        %compute the amplitude in the prestim period for real and pseudo
        split_real = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_real,:)),2));
        split_pseudo = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_pseudo,:)),2));
        %==fix me please               SD_split
        
        %for each channel, calculate the TTV pseudo and real relative to
        %their respective prestim high and low
        
        SD_pseudo = zeros(length(timefreq_data{q}.label),size(timefreq_data{q}.matrix,2),2);
        SD_real = SD_pseudo;
        for c = 1:length(timefreq_data{q}.label)
            splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
            SD_pseudo(c,:,1) = std(real(timefreq_data{q}.matrix(c,:,find(~splitindex))),[],3); %pseudo prestim low
            SD_pseudo(c,:,2) = std(real(timefreq_data{q}.matrix(c,:,find(splitindex))),[],3); %pseudo prestim high
            splitindex = split_real(c,:) > median(split_real(c,:));
            SD_real(c,:,1) = std(real(timefreq_data{q}.matrix(c,:,find(~splitindex))),[],3); %real prestim low
            SD_real(c,:,2) = std(real(timefreq_data{q}.matrix(c,:,find(splitindex))),[],3); %real prestim high
        end
        allmeas{q}.naddsd.amp.pseudo(:,:,:,i) = SD_pseudo;
        allmeas{q}.naddsd.amp.real(:,:,:,i) = SD_real;
        
        
        %normalize the TTV to its respective prestim
        allmeas{q}.naddttv.amp.pseudo(:,:,:,i) = (SD_pseudo(:,poststim_pseudo,:)-...
            mean(SD_pseudo(:,prestim_pseudo,:),2));
        allmeas{q}.naddttv.amp.real(:,:,:,i) = (SD_real(:,poststim_real,:)-...
            mean(SD_real(:,prestim_real,:),2));
        
        switch settings.units
            case 'prcchange'
                allmeas{q}.naddttv.amp.pseudo(:,:,:,i) = 100*allmeas{q}.naddttv.amp.pseudo(:,:,:,i)./...
                    (cat(3,mean(mean(SD_pseudo(:,prestim_pseudo,:),2),3),mean(mean(SD_pseudo(:,prestim_pseudo,:),2),3)));
                allmeas{q}.naddttv.amp.real(:,:,:,i) = 100*allmeas{q}.naddttv.amp.real(:,:,:,i)./...
                    (cat(3,mean(mean(SD_real(:,prestim_real,:),2),3),mean(mean(SD_real(:,prestim_real,:),2),3)));
            case 'zscore'
                allmeas{q}.naddttv.amp.pseudo(:,:,:,i) = zscore(allmeas{q}.naddttv.amp.pseudo(:,:,:,i),0,2);
                allmeas{q}.naddttv.amp.real(:,:,:,i) = zscore(allmeas{q}.naddttv.amp.real(:,:,:,i),0,2);
            case 'log'
                allmeas{q}.naddttv.amp.pseudo(:,:,:,i) = 10*log10(allmeas{q}.naddttv.amp.pseudo(:,:,:,i));
                allmeas{q}.naddttv.amp.real(:,:,:,i) = 10*log10(allmeas{q}.naddttv.amp.real(:,:,:,i));
        end
        %% Nonadditive TTV - phase-based
        split_real = squeeze(angle(timefreq_data{q}.matrix(:,prestim_real(end),:)));
        split_pseudo = squeeze(angle(timefreq_data{q}.matrix(:,prestim_pseudo(end),:)));
        
        %for each channel, calculate the TTV pseudo and real relative to
        %their respective prestim high and low
        
        SD_pseudo = zeros(length(timefreq_data{q}.label),size(timefreq_data{q}.matrix,2),2);
        SD_real = SD_pseudo;
        for c = 1:length(timefreq_data{q}.label)
            splitindex = cos(split_pseudo(c,:)) > median(cos(split_pseudo(c,:)));
            SD_pseudo(c,:,1) = std(real(timefreq_data{q}.matrix(c,:,find(~splitindex))),[],3); %pseudo prestim trough
            SD_pseudo(c,:,2) = std(real(timefreq_data{q}.matrix(c,:,find(splitindex))),[],3); %pseudo prestim peak
            
            splitindex = cos(split_real(c,:)) > median(cos(split_real(c,:)));
            SD_real(c,:,1) = std(real(timefreq_data{q}.matrix(c,:,find(~splitindex))),[],3); %real prestim trough
            SD_real(c,:,2) = std(real(timefreq_data{q}.matrix(c,:,find(splitindex))),[],3); %real prestim peak
        end
        allmeas{q}.naddsd.cos.pseudo(:,:,:,i) = SD_pseudo;
        allmeas{q}.naddsd.cos.real(:,:,:,i) = SD_real;
        
        
        %normalize the TTV to its respective prestim
        allmeas{q}.naddttv.cos.pseudo(:,:,:,i) = (SD_pseudo(:,poststim_pseudo,:)-...
            mean(SD_pseudo(:,prestim_pseudo,:),2));
        allmeas{q}.naddttv.cos.real(:,:,:,i) = (SD_real(:,poststim_real,:)-...
            mean(SD_real(:,prestim_real,:),2));
        
        switch settings.units
            case 'prcchange'
                allmeas{q}.naddttv.cos.pseudo(:,:,:,i) = 100*allmeas{q}.naddttv.cos.pseudo(:,:,:,i)./...
                    (cat(3,mean(mean(SD_pseudo(:,prestim_pseudo,:),2),3),mean(mean(SD_pseudo(:,prestim_pseudo,:),2),3)));
                allmeas{q}.naddttv.cos.real(:,:,:,i) = 100*allmeas{q}.naddttv.cos.real(:,:,:,i)./...
                    (cat(3,mean(mean(SD_real(:,prestim_real,:),2),3),mean(mean(SD_real(:,prestim_real,:),2),3)));
            case 'zscore'
                allmeas{q}.naddttv.cos.pseudo(:,:,:,i) = zscore(allmeas{q}.naddttv.cos.pseudo(:,:,:,i),0,2);
                allmeas{q}.naddttv.cos.real(:,:,:,i) = zscore(allmeas{q}.naddttv.cos.real(:,:,:,i),0,2);
            case 'log'
                allmeas{q}.naddttv.cos.pseudo(:,:,:,i) = 10*log10(allmeas{q}.naddttv.cos.pseudo(:,:,:,i));
                allmeas{q}.naddttv.cos.real(:,:,:,i) = 10*log10(allmeas{q}.naddttv.cos.real(:,:,:,i));
        end
        
        %% ERSP and ITC
        datacat = timefreq_data{q}.matrix;
        
        allmeas{q}.raw.ersp(:,:,i) = mean(abs(timefreq_data{q}.matrix),3);
        allmeas{q}.raw.itc(:,:,i) = abs(mean(datacat./abs(datacat),3));
        
        allmeas{q}.ersp.pseudo(:,:,i) = (mean(abs(datacat(:,poststim_pseudo,:)),3)...
            -mean(mean(abs(datacat(:,prestim_pseudo,:)),3),2));
        allmeas{q}.ersp.real(:,:,i) = (mean(abs(datacat(:,poststim_real,:)),3)...
            -mean(mean(abs(datacat(:,prestim_real,:)),3),2));
        
        switch settings.units
            case 'prcchange'
                allmeas{q}.ersp.pseudo(:,:,i) = 100*allmeas{q}.ersp.pseudo(:,:,i)./mean(mean(abs(datacat(:,prestim_pseudo,:)),3),2);
                allmeas{q}.ersp.real(:,:,i) = 100*allmeas{q}.ersp.real(:,:,i)./mean(mean(abs(datacat(:,prestim_real,:)),3),2);
            case 'zscore'
                allmeas{q}.ersp.pseudo(:,:,i) = zscore(allmeas{q}.ersp.pseudo(:,:,i),0,2);
                allmeas{q}.ersp.real(:,:,i) = zscore(allmeas{q}.ersp.real(:,:,i),0,2);
            case 'log'
                allmeas{q}.ersp.pseudo(:,:,i) = 10*log10(allmeas{q}.ersp.pseudo(:,:,i));
                allmeas{q}.ersp.real(:,:,i) = 10*log10(allmeas{q}.ersp.real(:,:,i));
        end
        
        allmeas{q}.itc.real(:,:,i) = (allmeas{q}.raw.itc(:,poststim_real,i)-mean(allmeas{q}.raw.itc(:,prestim_real,i),2));
        allmeas{q}.itc.pseudo(:,:,i) = (allmeas{q}.raw.itc(:,poststim_pseudo,i)-mean(allmeas{q}.raw.itc(:,prestim_pseudo,i),2))./...
            mean(allmeas{q}.raw.itc(:,prestim_pseudo,i),2);
        
        switch settings.units
            case 'prcchange'
                allmeas{q}.itc.real(:,:,i) = 100*allmeas{q}.itc.real(:,:,i)./mean(allmeas{q}.raw.itc(:,prestim_real,i),2);
                allmeas{q}.itc.pseudo(:,:,i) = 100*allmeas{q}.itc.pseudo(:,:,i)./mean(allmeas{q}.raw.itc(:,prestim_pseudo,i),2);
            case 'zscore'
                allmeas{q}.itc.real(:,:,i) = zscore(allmeas{q}.itc.real(:,:,i),0,2);
                allmeas{q}.itc.pseudo(:,:,i) = zscore(allmeas{q}.itc.pseudo(:,:,i),0,2);
            case 'log'
                allmeas{q}.itc.real(:,:,i) = 10*log10(allmeas{q}.itc.real(:,:,i));
                allmeas{q}.itc.pseudo(:,:,i) = 10*log10(allmeas{q}.itc.pseudo(:,:,i));
        end
        
        split_real = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_real,:)),2));
        split_pseudo = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_pseudo,:)),2));
        
        for c = 1:nbchan
            splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
            
            allmeas{q}.naddersp.amp.raw.pseudo(c,:,1,i) = mean(real(datacat(c,:,find(~splitindex))),3); %NOT ACTUALLY THE RAW ERSP - used for plotting later
            allmeas{q}.naddersp.amp.raw.pseudo(c,:,2,i) = mean(real(datacat(c,:,find(splitindex))),3);
            
            allmeas{q}.naddersp.amp.pseudo(c,:,1,i) = (mean(abs(datacat(c,poststim_pseudo,find(~splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),3),2)); %removed find(~splitindex) in the denominator because want to divide by a common baseline
            allmeas{q}.naddersp.amp.pseudo(c,:,2,i) = (mean(abs(datacat(c,poststim_pseudo,find(splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_pseudo,find(splitindex))),3),2));
            
            switch settings.units
                case 'prcchange'
                    allmeas{q}.naddersp.amp.pseudo(c,:,1,i) = 100*allmeas{q}.naddersp.amp.pseudo(c,:,1,i)./mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
                    allmeas{q}.naddersp.amp.pseudo(c,:,2,i) = 100*allmeas{q}.naddersp.amp.pseudo(c,:,2,i)./mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
                case 'zscore'
                    allmeas{q}.naddersp.amp.pseudo(c,:,1,i) = zscore(allmeas{q}.naddersp.amp.pseudo(c,:,1,i),0,2);
                    allmeas{q}.naddersp.amp.pseudo(c,:,2,i) = zscore(allmeas{q}.naddersp.amp.pseudo(c,:,2,i),0,2);
                case 'log'
                    allmeas{q}.naddersp.amp.pseudo(c,:,1,i) = 10*log10(allmeas{q}.naddersp.amp.pseudo(c,:,1,i));
                    allmeas{q}.naddersp.amp.pseudo(c,:,2,i) = 10*log10(allmeas{q}.naddersp.amp.pseudo(c,:,2,i))
            end
            
            splitindex = split_real(c,:) > median(split_real(c,:));
            
            allmeas{q}.naddersp.amp.raw.real(c,:,1,i) = mean(real(datacat(c,:,find(~splitindex))),3);
            allmeas{q}.naddersp.amp.raw.real(c,:,2,i) = mean(real(datacat(c,:,find(splitindex))),3);
            
            allmeas{q}.naddersp.amp.real(c,:,1,i) = (mean(abs(datacat(c,poststim_real,find(~splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_real,find(~splitindex))),3),2));
            allmeas{q}.naddersp.amp.real(c,:,2,i) = (mean(abs(datacat(c,poststim_real,find(splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_real,find(splitindex))),3),2));
            
            switch settings.units
                case 'prcchange'
                    allmeas{q}.naddersp.amp.real(c,:,1,i) = 100*allmeas{q}.naddersp.amp.real(c,:,1,i)./mean(mean(abs(datacat(c,prestim_real,:)),3),2);
                    allmeas{q}.naddersp.amp.real(c,:,2,i) = 100*allmeas{q}.naddersp.amp.real(c,:,2,i)./mean(mean(abs(datacat(c,prestim_real,:)),3),2);
                case 'zscore'
                    allmeas{q}.naddersp.amp.real(c,:,1,i) = zscore(allmeas{q}.naddersp.amp.real(c,:,1,i),0,2);
                    allmeas{q}.naddersp.amp.real(c,:,2,i) = zscore(allmeas{q}.naddersp.amp.real(c,:,2,i),0,2);
                case 'log'
                    allmeas{q}.naddersp.amp.real(c,:,1,i) = 10*log10(allmeas{q}.naddersp.amp.real(c,:,1,i));
                    allmeas{q}.naddersp.amp.real(c,:,2,i) = 10*log10(allmeas{q}.naddersp.amp.real(c,:,2,i));
            end
        end
        
        split_real = squeeze(angle(timefreq_data{q}.matrix(:,prestim_real(end),:)));
        split_pseudo = squeeze(angle(timefreq_data{q}.matrix(:,prestim_pseudo(end),:)));
        
        for c = 1:nbchan
            splitindex = cos(split_pseudo(c,:)) > median(cos(split_pseudo(c,:)));
            
            allmeas{q}.naddersp.cos.raw.pseudo(c,:,1,i) = mean(real(datacat(c,:,find(~splitindex))),3);
            allmeas{q}.naddersp.cos.raw.pseudo(c,:,2,i) = mean(real(datacat(c,:,find(splitindex))),3);
            
            allmeas{q}.naddersp.cos.pseudo(c,:,1,i) = (mean(abs(datacat(c,poststim_pseudo,find(~splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),3),2));
            allmeas{q}.naddersp.cos.pseudo(c,:,2,i) = (mean(abs(datacat(c,poststim_pseudo,find(splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_pseudo,find(splitindex))),3),2));
            
            switch settings.units
                case 'prcchange'
                    allmeas{q}.naddersp.cos.pseudo(c,:,1,i) = 100*allmeas{q}.naddersp.cos.pseudo(c,:,1,i)./mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
                    allmeas{q}.naddersp.cos.pseudo(c,:,2,i) = 100*allmeas{q}.naddersp.cos.pseudo(c,:,2,i)./mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
                case 'zscore'
                    allmeas{q}.naddersp.cos.pseudo(c,:,1,i) = zscore(allmeas{q}.naddersp.cos.pseudo(c,:,1,i),0,2);
                    allmeas{q}.naddersp.cos.pseudo(c,:,2,i) = zscore(allmeas{q}.naddersp.cos.pseudo(c,:,2,i),0,2);
                case 'log'
                    allmeas{q}.naddersp.cos.real(c,:,1,i) = 10*log10(allmeas{q}.naddersp.cos.pseudo(c,:,1,i));
                    allmeas{q}.naddersp.cos.pseudo(c,:,2,i) = 10*log10(allmeas{q}.naddersp.cos.pseudo(c,:,2,i));
            end
            
            splitindex = cos(split_real(c,:)) > median(cos(split_real(c,:)));
            
            allmeas{q}.naddersp.cos.raw.real(c,:,1,i) = mean(real(datacat(c,:,find(~splitindex))),3);
            allmeas{q}.naddersp.cos.raw.real(c,:,2,i) = mean(real(datacat(c,:,find(splitindex))),3);
            
            allmeas{q}.naddersp.cos.real(c,:,1,i) = (mean(abs(datacat(c,poststim_real,find(~splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_real,find(~splitindex))),3),2));
            allmeas{q}.naddersp.cos.real(c,:,2,i) = (mean(abs(datacat(c,poststim_real,find(splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_real,find(splitindex))),3),2));
            
            switch settings.units
                case 'prcchange'
                    allmeas{q}.naddersp.cos.real(c,:,1,i) = 100*allmeas{q}.naddersp.cos.real(c,:,1,i)./mean(mean(abs(datacat(c,prestim_real,:)),3),2);
                    allmeas{q}.naddersp.cos.real(c,:,2,i) = 100*allmeas{q}.naddersp.cos.real(c,:,2,i)./mean(mean(abs(datacat(c,prestim_real,:)),3),2);
                case 'zscore'
                    allmeas{q}.naddersp.cos.real(c,:,1,i) = zscore(allmeas{q}.naddersp.cos.real(c,:,1,i),0,2);
                    allmeas{q}.naddersp.cos.real(c,:,2,i) = zscore(allmeas{q}.naddersp.cos.real(c,:,2,i),0,2);
                case 'log'
                    allmeas{q}.naddersp.cos.real(c,:,1,i) = 10*log10(allmeas{q}.naddersp.cos.real(c,:,1,i));
                    allmeas{q}.naddersp.cos.real(c,:,2,i) = 10*log10(allmeas{q}.naddersp.cos.real(c,:,2,i))
            end
        end
        %% TTV of ERSP
        
        allmeas{q}.raw.ttversp(:,:,i) = std(abs(datacat),[],3);
        allmeas{q}.ttversp.pseudo(:,:,i) = (std(abs(datacat(:,poststim_pseudo,:)),[],3)...
            -mean(std(abs(datacat(:,prestim_pseudo,:)),[],3),2));
        allmeas{q}.ttversp.real(:,:,i) = (std(abs(datacat(:,poststim_real,:)),[],3)...
            -mean(std(abs(datacat(:,prestim_real,:)),[],3),2));
        
        switch settings.units
            case 'prcchange'
                allmeas{q}.ttversp.real(:,:,i) = 100*allmeas{q}.ttversp.real(:,:,i)./mean(allmeas{q}.raw.ttversp(:,prestim_real,i),2);
                allmeas{q}.ttversp.pseudo(:,:,i) = 100*allmeas{q}.ttversp.pseudo(:,:,i)./mean(allmeas{q}.raw.ttversp(:,prestim_pseudo,i),2);
            case 'zscore'
                allmeas{q}.ttversp.real(:,:,i) = zscore(allmeas{q}.ttversp.real(:,:,i),0,2);
                allmeas{q}.ttversp.pseudo(:,:,i) = zscore(allmeas{q}.ttversp.pseudo(:,:,i),0,2);
            case 'log'
                allmeas{q}.ttversp.real(:,:,i) = 10*log10(allmeas{q}.ttversp.real(:,:,i));
                allmeas{q}.ttversp.pseudo(:,:,i) = 10*log10(allmeas{q}.ttversp.pseudo(:,:,i));
        end
        
        %% ERP nonadditivity
        allmeas{q}.erp(:,:,i) = mean(real(timefreq_data{q}.matrix),3);
        
        split_real = squeeze(mean(real(timefreq_data{q}.matrix(:,prestim_real,:)),2));
        split_pseudo = squeeze(mean(real(timefreq_data{q}.matrix(:,prestim_pseudo,:)),2));
        
        %datacat = timefreq_data{q}.matrix;
        
        for c = 1:nbchan
            splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
            
            allmeas{q}.nadderp.raw.pseudo(c,:,1,i) = mean(real(datacat(c,:,find(~splitindex))),3); 
            allmeas{q}.nadderp.raw.pseudo(c,:,2,i) = mean(real(datacat(c,:,find(splitindex))),3);
            
            
            allmeas{q}.nadderp.pseudo(c,:,1,i) = (mean(real(datacat(c,poststim_pseudo,find(~splitindex))),3)...
                -mean(mean(real(datacat(c,prestim_pseudo,find(~splitindex))),3),2));
            allmeas{q}.nadderp.pseudo(c,:,2,i) = (mean(real(datacat(c,poststim_pseudo,find(splitindex))),3)...
                -mean(mean(real(datacat(c,prestim_pseudo,find(splitindex))),3),2));
            
%             switch settings.units
%                 case 'prcchange'
%                     allmeas{q}.nadderp.pseudo(c,:,1,i) = 100*allmeas{q}.nadderp.pseudo(c,:,1,i)./mean(mean(real(datacat(c,prestim_pseudo,:)),3),2);
%                     allmeas{q}.nadderp.pseudo(c,:,2,i) = 100*allmeas{q}.nadderp.pseudo(c,:,2,i)./mean(mean(real(datacat(c,prestim_pseudo,:)),3),2);
%                 case 'zscore'
%                     allmeas{q}.nadderp.pseudo(c,:,1,i) = zscore(allmeas{q}.nadderp.pseudo(c,:,1,i),0,2);
%                     allmeas{q}.nadderp.pseudo(c,:,2,i) = zscore(allmeas{q}.nadderp.pseudo(c,:,2,i),0,2);
%                 case 'log'
%                     allmeas{q}.nadderp.pseudo(c,:,1,i) = 10*log10(allmeas{q}.nadderp.pseudo(c,:,1,i));
%                     allmeas{q}.nadderp.pseudo(c,:,2,i) = 10*log10(allmeas{q}.nadderp.pseudo(c,:,2,i));
%             end
            
            splitindex = split_real(c,:) > median(split_real(c,:));
            
            allmeas{q}.nadderp.raw.real(c,:,1,i) = mean(real(datacat(c,:,find(~splitindex))),3); %NOT ACTUALLY THE RAW ERSP - used for plotting later
            allmeas{q}.nadderp.raw.real(c,:,2,i) = mean(real(datacat(c,:,find(splitindex))),3);
            
            
            allmeas{q}.nadderp.real(c,:,1,i) = (mean(real(datacat(c,poststim_real,find(~splitindex))),3)...
                -mean(mean(real(datacat(c,prestim_real,find(~splitindex))),3),2));
            allmeas{q}.nadderp.real(c,:,2,i) = (mean(real(datacat(c,poststim_real,find(splitindex))),3)...
                -mean(mean(real(datacat(c,prestim_real,find(splitindex))),3),2));
                    
%             switch settings.units
%                 case 'prcchange'
%                     allmeas{q}.nadderp.real(c,:,1,i) = 100*allmeas{q}.nadderp.real(c,:,1,i)./mean(mean(real(datacat(c,prestim_real,:)),3),2);
%                     allmeas{q}.nadderp.real(c,:,2,i) = 100*allmeas{q}.nadderp.real(c,:,2,i)./mean(mean(real(datacat(c,prestim_real,:)),3),2);
%                 case 'zscore'
%                     allmeas{q}.nadderp.real(c,:,1,i) = zscore(allmeas{q}.nadderp.real(c,:,1,i),0,2);
%                     allmeas{q}.nadderp.real(c,:,2,i) = zscore(allmeas{q}.nadderp.real(c,:,2,i),0,2);
%                 case 'log'
%                     allmeas{q}.nadderp.real(c,:,1,i) = 10*log10(allmeas{q}.nadderp.real(c,:,1,i));
%                     allmeas{q}.nadderp.real(c,:,2,i) = 10*log10(allmeas{q}.nadderp.real(c,:,2,i))
%             end
        end
        
        
    end
end

if ~strcmpi(settings.tfparams.method,'hilbert')
    timefreq_data = parload(files(1).name,'timefreq_data');
    for c = 2:length(timefreq_data)
       parents(c) = timefreq_data{c}.parent;
    end
    newmeas = cell(1,length(settings.tfparams.fbands));
    for c = 2:length(settings.tfparams.fbands)
        children{c} = find(parents == c);
        fields = fieldnames_recurse(allmeas{c});
        fields = cell_unpack(fields);
       for cc = 1:length(fields)
           tmp = [];
           for ccc = 1:length(children{c})
               dimn = length(size(getfield_nest(allmeas{children{c}(ccc)},fields{cc})));
               tmp = cat(dimn+1,tmp,getfield_nest(allmeas{children{c}(ccc)},(fields{cc})));
           end
           newmeas{c} = assignfield_nest(newmeas{c},fields{cc},nanmean(tmp,dimn+1));
       end
    end
    
    dimn = [];
    children{1} = 2:length(timefreq_data);
    fields = fieldnames_recurse(allmeas{1});
    fields = cell_unpack(fields);

       for cc = 1:length(fields)
           tmp = [];
           for ccc = 1:length(children{1})
               dimn = length(size(getfield_nest(allmeas{children{1}(ccc)},fields{cc})));
               tmp = cat(dimn+1,tmp,getfield_nest(allmeas{children{1}(ccc)},(fields{cc})));
           end
           newmeas{1} = assignfield_nest(newmeas{1},fields{cc},nanmean(tmp,dimn+1));
       end
    newmeas{1}.erp = allmeas{1}.erp;
    newmeas{1}.nadderp = allmeas{1}.nadderp;
    newmeas{1}.ttv = allmeas{1}.ttv;
    
    settings.nfreqs = length(settings.tfparams.fbands);
    numbands = settings.nfreqs;
    
    allmeas = newmeas;
    newmeas = [];
end


%% Calculating indices
for q = 1:numbands
    allmeas{q}.ttvindex = squeeze(trapz(allmeas{q}.ttv.real(:,aucindex,:),2));% - squeeze(trapz(allmeas{q}.ttv.pseudo(:,aucindex,:),2));
    allmeas{q}.erspindex = squeeze(trapz(allmeas{q}.ersp.real(:,aucindex,:),2));% - squeeze(trapz(allmeas{q}.ersp.pseudo(:,aucindex,:),2));
    allmeas{q}.itcindex = squeeze(trapz(allmeas{q}.itc.real(:,aucindex,:),2));% - squeeze(trapz(allmeas{q}.itc.pseudo(:,aucindex,:),2));
    allmeas{q}.ttverspindex = squeeze(trapz(allmeas{q}.ttversp.real(:,aucindex,:),2));% - squeeze(trapz(allmeas{q}.ttversp.pseudo(:,aucindex,:),2));
    
    % abs((prestim high real - pseudo)) - abs((prestim low real - pseudo))
    allmeas{q}.nattvindex.amp = abs((squeeze(trapz(allmeas{q}.naddttv.amp.real(:,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.naddttv.amp.pseudo(:,aucindex,2,:),2)))) - ...
        abs((squeeze(trapz(allmeas{q}.naddttv.amp.real(:,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.naddttv.amp.pseudo(:,aucindex,1,:),2))));
    allmeas{q}.naerspindex.amp = abs((squeeze(trapz(allmeas{q}.naddersp.amp.real(:,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.naddersp.amp.pseudo(:,aucindex,2,:),2)))) - ...
        abs((squeeze(trapz(allmeas{q}.naddersp.amp.real(:,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.naddersp.amp.pseudo(:,aucindex,1,:),2))));
    
    allmeas{q}.nattvindex.cos = abs((squeeze(trapz(allmeas{q}.naddttv.cos.real(:,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.naddttv.cos.pseudo(:,aucindex,2,:),2)))) - ...
        abs((squeeze(trapz(allmeas{q}.naddttv.cos.real(:,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.naddttv.cos.pseudo(:,aucindex,1,:),2))));
    allmeas{q}.naerspindex.cos = abs((squeeze(trapz(allmeas{q}.naddersp.cos.real(:,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.naddersp.cos.pseudo(:,aucindex,2,:),2)))) - ...
        abs((squeeze(trapz(allmeas{q}.naddersp.cos.real(:,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.naddersp.cos.pseudo(:,aucindex,1,:),2))));
    
    allmeas{q}.naerpindex = abs((squeeze(trapz(allmeas{q}.nadderp.real(:,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.nadderp.pseudo(:,aucindex,2,:),2)))) - ...
        abs((squeeze(trapz(allmeas{q}.nadderp.real(:,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.nadderp.pseudo(:,aucindex,1,:),2))));
end

settings.pseudo.prestim = prestim_pseudo;
settings.real.prestim = prestim_real;
settings.pseudo.poststim = poststim_pseudo;
settings.real.poststim = poststim_real;

%% Saving and cleaning up
save([settings.outputdir '/' settings.datasetname '_calc.mat'],'allmeas','-v7.3')

delete(gcp('nocreate'))
end