function TTV_ERSP_calc_func(settings)

numbands = length(settings.tfparams.fbandnames);

fbands = settings.tfparams.fbandnames;

aucindex = settings.aucindex;

if isempty(gcp('nocreate'))
    parpool(numbands)
end

files = dir('*freq_filtered*');

allmeas = cell(1,length(fbands));
for c = 1:7
   allmeas{c} = struct; 
end

prestim_pseudo = settings.pseudo.prestim;
prestim_real = settings.real.prestim;
poststim_pseudo = settings.pseudo.poststim;
poststim_real = settings.real.poststim;

for i = 1:length(files)
    disp(['Processing subject ' num2str(i)])
    timefreq_data = parload(files(i).name,'timefreq_data');
    nbchan = size(timefreq_data{1}.trial{1},1);
    parfor q = 1:numbands
        timefreq_data{q}.matrix = cat(3,timefreq_data{q}.trial{:});
        
        SD_all = zeros(length(timefreq_data{q}.label),size(timefreq_data{q}.matrix,2));
        for c = 1:length(timefreq_data{q}.label)
            SD_all(c,:) = std(real(timefreq_data{q}.matrix(c,:,:)),[],3);
            %TTV_real(c,:) = std(read(timefreq_data{q}.matrix(c,191:380,:)),[],3);
        end
        allmeas{q}.raw.sd(:,:,i) = SD_all;
        
        % normalize by 50ms prestim ifor both real and pseudotrial
        allmeas{q}.ttv.pseudo(:,:,i) = 100*(SD_all(:,poststim_pseudo)-...
            mean(SD_all(:,prestim_pseudo),2))./mean(SD_all(:,prestim_pseudo),2);
        allmeas{q}.ttv.real(:,:,i) = 100*(SD_all(:,poststim_real)-...
            mean(SD_all(:,prestim_real),2))./mean(SD_all(:,prestim_real),2);
        
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
        allmeas{q}.naddttv.amp.pseudo(:,:,:,i) = 100*(SD_pseudo(:,poststim_pseudo,:)-...
            mean(SD_pseudo(:,prestim_pseudo,:),2))./(cat(3,mean(mean(SD_pseudo(:,prestim_pseudo,:),2),3),...
            mean(mean(SD_pseudo(:,prestim_pseudo,:),2),3)));
        allmeas{q}.naddttv.amp.real(:,:,:,i) = 100*(SD_real(:,poststim_real,:)-...
            mean(SD_real(:,prestim_real,:),2))./(cat(3,mean(mean(SD_real(:,prestim_real,:),2),3),...
            mean(mean(SD_real(:,prestim_real,:),2),3)));
        
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
        allmeas{q}.naddttv.cos.pseudo(:,:,:,i) = 100*(SD_pseudo(:,poststim_pseudo,:)-...
            mean(SD_pseudo(:,prestim_pseudo,:),2))./(cat(3,mean(mean(SD_pseudo(:,prestim_pseudo,:),2),3),...
            mean(mean(SD_pseudo(:,prestim_pseudo,:),2),3)));
        allmeas{q}.naddttv.cos.real(:,:,:,i) = 100*(SD_real(:,poststim_real,:)-...
            mean(SD_real(:,prestim_real,:),2))./(cat(3,mean(mean(SD_real(:,prestim_real,:),2),3),...
            mean(mean(SD_real(:,prestim_real,:),2),3)));
        
        %% ERSP and ITC
        datacat = timefreq_data{q}.matrix;

        allmeas{q}.raw.ersp(:,:,i) = mean(abs(timefreq_data{q}.matrix),3);
        allmeas{q}.raw.itc(:,:,i) = abs(mean(datacat./abs(datacat),3));
                
        allmeas{q}.ersp.pseudo(:,:,i) = 100*(mean(abs(datacat(:,poststim_pseudo,:)),3)...
            -mean(mean(abs(datacat(:,prestim_pseudo,:)),3),2))./mean(mean(abs(datacat(:,prestim_pseudo,:)),3),2);
        allmeas{q}.ersp.real(:,:,i) = 100*(mean(abs(datacat(:,poststim_real,:)),3)...
            -mean(mean(abs(datacat(:,prestim_real,:)),3),2))./mean(mean(abs(datacat(:,prestim_real,:)),3),2);
        
        allmeas{q}.itc.real(:,:,i) = 100*(allmeas{q}.raw.itc(:,poststim_real,i)-mean(allmeas{q}.raw.itc(:,prestim_real,i),2))./...
            mean(allmeas{q}.raw.itc(:,prestim_real,i),2);
        allmeas{q}.itc.pseudo(:,:,i) = 100*(allmeas{q}.raw.itc(:,poststim_pseudo,i)-mean(allmeas{q}.raw.itc(:,prestim_pseudo,i),2))./...
            mean(allmeas{q}.raw.itc(:,prestim_pseudo,i),2);
        
        split_real = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_real,:)),2));
        split_pseudo = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_pseudo,:)),2));
        
        for c = 1:nbchan
            splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
            
            allmeas{q}.naddersp.amp.pseudo(c,:,1,i) = 100*(mean(abs(datacat(c,poststim_pseudo,find(~splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),3),2))./...
                mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);  %removed find(~splitindex) in the denominator because want to divide by a common baseline
            allmeas{q}.naddersp.amp.pseudo(c,:,2,i) = 100*(mean(abs(datacat(c,poststim_pseudo,find(splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_pseudo,find(splitindex))),3),2))./...
                mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
            
            splitindex = split_real(c,:) > median(split_real(c,:));
            
            allmeas{q}.naddersp.amp.real(c,:,1,i) = 100*(mean(abs(datacat(c,poststim_real,find(~splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_real,find(~splitindex))),3),2))./...
                mean(mean(abs(datacat(c,prestim_real,:)),3),2);
            allmeas{q}.naddersp.amp.real(c,:,2,i) = 100*(mean(abs(datacat(c,poststim_real,find(splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_real,find(splitindex))),3),2))./...
                mean(mean(abs(datacat(c,prestim_real,:)),3),2);
        end
        
        split_real = squeeze(angle(timefreq_data{q}.matrix(:,800,:)));
        split_pseudo = squeeze(angle(timefreq_data{q}.matrix(:,600,:)));
        
        for c = 1:nbchan
            splitindex = cos(split_pseudo(c,:)) > median(cos(split_pseudo(c,:)));
            
            allmeas{q}.naddersp.cos.pseudo(c,:,1,i) = 100*(mean(abs(datacat(c,poststim_pseudo,find(~splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),3),2))./...
                mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
            allmeas{q}.naddersp.cos.pseudo(c,:,2,i) = 100*(mean(abs(datacat(c,poststim_pseudo,find(splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_pseudo,find(splitindex))),3),2))./...
                mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
            
            splitindex = cos(split_real(c,:)) > median(cos(split_real(c,:)));
            
            allmeas{q}.naddersp.cos.real(c,:,1,i) = 100*(mean(abs(datacat(c,poststim_real,find(~splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_real,find(~splitindex))),3),2))./...
                mean(mean(abs(datacat(c,prestim_real,:)),3),2);
            allmeas{q}.naddersp.cos.real(c,:,2,i) = 100*(mean(abs(datacat(c,poststim_real,find(splitindex))),3)...
                -mean(mean(abs(datacat(c,prestim_real,find(splitindex))),3),2))./...
                mean(mean(abs(datacat(c,prestim_real,:)),3),2);
        end
    end
end


%% Calculating indices
for q = 1:numbands
    allmeas{q}.ttvindex = squeeze(trapz(allmeas{q}.ttv.real(:,aucindex,:),2)) - squeeze(trapz(allmeas{q}.ttv.pseudo(:,aucindex,:),2));
    allmeas{q}.erspindex = squeeze(trapz(allmeas{q}.ersp.real(:,aucindex,:),2)) - squeeze(trapz(allmeas{q}.ersp.pseudo(:,aucindex,:),2));
    allmeas{q}.itcindex = squeeze(trapz(allmeas{q}.itc.real(:,aucindex,:),2)) - squeeze(trapz(allmeas{q}.itc.pseudo(:,aucindex,:),2));
    
    % (prestim high real - pseudo) - (prestim low real - pseudo)
    allmeas{q}.nattvindex.amp = (squeeze(trapz(allmeas{q}.naddttv.amp.real(:,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.naddttv.amp.pseudo(:,aucindex,2,:),2))) - ...
        (squeeze(trapz(allmeas{q}.naddttv.amp.real(:,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.naddttv.amp.pseudo(:,aucindex,1,:),2)));
    allmeas{q}.naerspindex.amp = (squeeze(trapz(allmeas{q}.naddersp.amp.real(:,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.naddersp.amp.pseudo(:,aucindex,2,:),2))) - ...
        (squeeze(trapz(allmeas{q}.naddersp.amp.real(:,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.naddersp.amp.pseudo(:,aucindex,1,:),2)));
    
    allmeas{q}.nattvindex.cos = (squeeze(trapz(allmeas{q}.naddttv.cos.real(:,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.naddttv.cos.pseudo(:,aucindex,2,:),2))) - ...
        (squeeze(trapz(allmeas{q}.naddttv.cos.real(:,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.naddttv.cos.pseudo(:,aucindex,1,:),2)));
    allmeas{q}.naerspindex.cos = (squeeze(trapz(allmeas{q}.naddersp.cos.real(:,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.naddersp.cos.pseudo(:,aucindex,2,:),2))) - ...
        (squeeze(trapz(allmeas{q}.naddersp.cos.real(:,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.naddersp.cos.pseudo(:,aucindex,1,:),2)));
end


%% Saving and cleaning up
save([settings.outputdir '/' settings.datasetname '_calc.mat'],'allmeas','-v7.3')

delete(gcp('nocreate'))
end