function [stats,indvar_record,depvar_record] = TrialCorr_Ivar(trange,frange,electrode,indvar_name,depvar_name,varargin)
%logistic regression of prestimulus phase on accuracy

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/trialFlags_ivar_2.mat')
trialFlags = rmfield(trialFlags,'sub214');
trialFlags = rmfield(trialFlags,'sub224');
if ~EasyParse(varargin,'LongEpochs','true')
    directory_name = '/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/raw_eeg/self/ICA_componentsrejected/All_eegcaps/';
else
    directory_name = '/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/raw_eeg/self/preprocessed_onlyContinuous/SelfICAMARA_LongEpochs';
end

files = dir(fullfile(directory_name, '*.mat'));

fileindex = find(~[files.isdir]);

filenames = extractfield(files,'name');

currentData = [];

fields = fieldnames(trialFlags);

currentData = [];

%meandataindex = find(meandata.subid == str2num(subID),1);

%if meandata.fnirscap(meandataindex)
elecindex.frontalf = [13 14 19 20];
elecindex.parietalf = [40 51 45 46];
%else
elecindex.frontaln = [37 38 9 10];
elecindex.parietaln = [52 24 56 57];

elecnames = cell(1,63);
%elecnames(elecindex.frontalf) = {'FFC1h','FFC2h','Fc1','Fc2'};
%elecnames(elecindex.parietalf) = {'Cpz','Pz','CPP1h','CPP2h'};

elecnames(elecindex.frontalf) = {'F1','F2','Fc1','Fc2'};
elecnames(elecindex.parietalf) = {'Cpz','Pz','P1','P2'};

elecnames(elecindex.frontaln) = {'F1','F2','Fc1','Fc2'};
elecnames(elecindex.parietaln) = {'Cpz','Pz','P1','P2'};

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/allself.mat');

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/meandata_2.mat');

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','TTV');

if ~EasyParse(varargin,'LongEpochs','true')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TFparams.mat');
else
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TFparams_long.mat');
end

orig_trange = trange;
orig_frange = frange;

[~,t1] = min(abs(timesTF-trange(1)));
[~,t2] = min(abs(timesTF-trange(2)));
trange = t1:t2;
disp(['Taking time points ' num2str(t1) ' to ' num2str(t2) '...'])

[~,f1] = min(abs(freqsTF-frange(1)));
[~,f2] = min(abs(freqsTF-frange(2)));
frange = f1:f2;
disp(['Taking frequencies ' num2str(f1) ' to ' num2str(f2) '...'])

indvar2 = [];
depvar2 = [];
if ~EasyParse(varargin,'Subjectbased','false')
    indvar_record = struct;
    depvar_record = struct;
end

if strcmpi(indvar_name,'difmean') || strcmpi(indvar_name,'amp_win') || strcmpi(indvar_name,'amplitude') || strcmpi(indvar_name,'amp_slope')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','poststim');
    %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTV_response.mat','response');
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/meandata_2.mat');
    trange = ceil(orig_trange/2);
    trange_ttv = trange(1):trange(2);
end

if strcmpi(indvar_name,'amp_split') || strcmpi(indvar_name,'amp_split_win')
    if (ischar(orig_frange) || isstring(orig_frange)) && ~strcmpi(orig_frange,'broadband')
           load(['/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTV' orig_frange '.mat'],'prestim');
    else
   load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','prestim');
    end
   prestim2 = prestim;
   clear prestim
   warning('Using pseudotrial')
   %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat','prestim2');
   trange_ttv = ceil(orig_trange/2);
   trange_ttv = trange_ttv(1):trange_ttv(2);
end


for q = 1:length(fields)
    
    subID = extractAfter(fields{q},'sub');
    
    if CheckInput(varargin,'Conditions')
        [~,~,flatindices] = Ivar_findConds_2(EasyParse(varargin,'Conditions'),trialFlags.(fields{q}));
        flatindices = flatindices{1};
    else
        flatindices = 1:400;
    end
    
    if strcmpi(indvar_name,'difmean') || strcmpi(indvar_name,'amp_win') || strcmpi(indvar_name,'amp_slope')
        fnirscap = meandata.fnirscap(find(meandata.subid == str2num(subID)));
        ttv_electrodes = zeros(1,length(electrode));
        
        ttv_elecindex.fnirs = [13 14 19 20 40 51 45 46];
        
        ttv_elecindex.normal = [37 38 9 10 52 24 56 57];
        
        ttv_elecnames.fnirs = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};
        
        ttv_elecnames.normal = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};
        
        if fnirscap
            poststim.(['sub' subID]) = squeeze(mean(poststim.(['sub' subID])(ttv_elecindex.fnirs(find(contains(ttv_elecnames.fnirs,electrode))),:,:),1));
        else
            poststim.(['sub' subID]) = squeeze(mean(poststim.(['sub' subID])(ttv_elecindex.normal(find(contains(ttv_elecnames.normal,electrode))),:,:),1));
        end
    end
    
    if strcmpi(indvar_name,'amp_split') || strcmpi(indvar_name,'amp_split_win')
        fnirscap = meandata.fnirscap(find(meandata.subid == str2num(subID)));
        ttv_electrodes = zeros(1,length(electrode));
        
        ttv_elecindex.fnirs = [13 14 19 20 40 51 45 46];
        
        ttv_elecindex.normal = [37 38 9 10 52 24 56 57];
        
        ttv_elecnames.fnirs = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};
        
        ttv_elecnames.normal = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};
        if iscell(electrode)
            if fnirscap
                prestim2.(['sub' subID]) = squeeze(mean(prestim2.(['sub' subID])(ttv_elecindex.fnirs(find(contains(ttv_elecnames.fnirs,electrode))),:,:),1));
            else
                prestim2.(['sub' subID]) = squeeze(mean(prestim2.(['sub' subID])(ttv_elecindex.normal(find(contains(ttv_elecnames.normal,electrode))),:,:),1));
            end
        else
            prestim2.(['sub' subID]) = squeeze(mean(prestim2.(['sub' subID]),1));
        end
    end
    
    newindex = contains(filenames,extractAfter(fields{q},'sub'));
    
    newindex = fileindex(find(newindex));
    
    if ~isempty(newindex)
        if strcmpi(indvar_name,'ersp') || strcmpi(indvar_name,'phase') || strcmpi(depvar_name,'ersp') || strcmpi(depvar_name,'phase') || strcmpi(indvar_name,'ersp_zscore') || strcmpi(indvar_name,'raw_ersp') || strcmpi(indvar_name,'onset_phase') || strcmpi(depvar_name,'ersp_timecourse') || strcmpi(indvar_name,'erspindex')
            
            for i = 1:length(newindex)
                
                filename = files(newindex(i)).name;
                [PATH, NAME, EXT] = fileparts(filename);
                
                NAME = [NAME, EXT];
                
                subID = extractBetween(filename,'tfdata_','_electrode');
                subID = subID{1};
                
                currElectrode = extractBetween(filename,'electrode','.mat');
                currElectrode = currElectrode{1};
                
                if any(contains(electrode,elecnames{str2num(currElectrode)}))
                    
                    load([directory_name '/' filename])
                    
                    try
                        %[~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(['sub' subID]));
                        currentData = cat(4,currentData,tfdata(frange,:,flatindices));
                        warning('When using multiple electrodes make sure to check your results')
                    catch
                        disp(lasterror)
                    end
                end
            end
            currentData = squeeze(currentData);
            
%             currentData = mean(mean(currentData,1),2);
%             currentData = squeeze(currentData);
%             currentData = reshape(currentData,400,length(electrode));
%             currentData = mean(currentData,2);
        end
        
        rmtrials = [];
        depvar = [];
        indvar = [];
        
        if strcmpi(depvar_name,'accuracy') || strcmpi(depvar_name,'acc_ttest')
            for c = 1:length(trialFlags.(['sub' subID]))
                if EasyParse(varargin,'AccFilter','on')
                    if trialFlags.(['sub' subID]){c}(4) == '3'
                        depvar(c) = 1;
                        %                 else
                        %                     depvar(c) = 0;
                        %                 end
                    elseif trialFlags.(['sub' subID]){c}(4) == '2'
                        depvar(c) = 0;
                    else
                        depvar(c) = NaN;
                        rmtrials = [rmtrials c];
                    end
                else
                    if trialFlags.(['sub' subID]){c}(4) == '3'
                        depvar(c) = 1;
                    else
                        depvar(c) = 0;
                    end
                end
            end
            depvar = depvar';
            depvar = depvar(flatindices);
            rmtrials = find(ismember(flatindices,rmtrials));
        elseif strcmpi(depvar_name,'RT')
            startPoint = find(allself.subid == str2num(subID),1);
            depvar = allself.reac_time(startPoint:startPoint+399);
            depvar = depvar(flatindices);
            rmtrials = find(depvar == -1)';
        elseif strcmpi(depvar_name,'ersp')
            depvar = mean(mean(bsxfun(@minus,10*log10(abs(currentData(:,trange,:)).^2),mean(10*log10(abs(currentData(:,1:60,:)).^2),2)),1),2); 
            depvar = squeeze(depvar);
        elseif strcmpi(depvar_name,'ersp_timecourse')
            depvar = mean(mean(bsxfun(@minus,10*log10(abs(currentData(:,trange,:)).^2),mean(10*log10(abs(currentData(:,1:60,:)).^2),2)),1),3); 
            depvar = squeeze(depvar);
        elseif strcmpi(depvar_name,'ersp_timecourse2')
            depvar = mean(bsxfun(@minus,10*log10(abs(currentData(:,trange,:)).^2),mean(10*log10(abs(currentData(:,1:60,:)).^2),2)),1);
            depvar = squeeze(depvar);
        end

        
        if strcmpi(indvar_name,'phase')
            %sinangle = sin(angle(currentData));
            %cosangle = cos(angle(currentData));
            indvar = circ_mean(circ_mean(angle(currentData),[],2),[],1);
            indvar = squeeze(indvar);
        elseif strcmpi(indvar_name,'ersp')
            for c = 1:size(currentData,4)
            indvar{c} = mean(mean(bsxfun(@minus,10*log10(abs(currentData(:,trange,:,c)).^2),mean(10*log10(abs(currentData(:,1:60,:,c)).^2),2)),1),2); 
            end
            indvar = squeeze(mean(cell2mat(indvar)));
            %indvar = 10*log10(abs(currentData).^2); %not baseline corrected
%         elseif strcmpi(indvar_name,'ttv')
%             if orig_trange(1) < 0
%                indvar = mean(TTV.prestim.(['sub' subID])(:,1000+(orig_trange(1)/2):(1000+origtrange(2)/2)),1);
%                indvar = squeeze(indvar);
%                indvar = indvar';
%             end
        elseif strcmpi(indvar_name,'erspindex')
            for c = 1:size(currentData,4)
            indvar{c} = mean(bsxfun(@minus,10*log10(abs(currentData(:,trange,:,c)).^2),mean(10*log10(abs(currentData(:,1:60,:,c)).^2),2)),1); 
            indvar{c} = trapz(indvar{c});
            end
            indvar = cell2mat(indvar);
            indvar = squeeze(mean(indvar));
        elseif strcmpi(indvar_name,'raw_ersp')
            indvar = mean(mean(10*log10(abs(currentData(:,52:60,:).^2))));
            indvar = squeeze(indvar);
            warning('Hard-coded time')
            
        elseif strcmpi(indvar_name,'onset_phase')
            indvar = circ_mean(circ_mean(angle(currentData(:,60,:)),[],2),[],1);
            indvar = squeeze(indvar);
            warning('Hard-coded time')
        elseif strcmpi(indvar_name,'difmean')
            indvar = poststim.(['sub' subID])(trange_ttv,:) - mean(poststim.(['sub' subID])(trange_ttv,:),2);
            indvar = mean(indvar,1);
            indvar = indvar';
            indvar = double(indvar);
        elseif strcmpi(indvar_name,'ersp_zscore')
            indvar = 10*log10(abs(currentData).^2); %not baseline corrected
            indvar = abs(zscore(indvar));
        elseif strcmpi(indvar_name,'amp_win')
            tmp = hilbert(poststim.(['sub' subID]));
            %submeanRT = meandata.MeanRT(find(meandata.subid == str2num(subID)));
%             startpoint = find(allself.subid == str2num(subID),1);
            %trange_ttv = ceil((submeanRT:(submeanRT+400))/2);
            tmp = abs(tmp(trange_ttv,:));
%             tmp = abs(tmp(:,:));
            %tmp = poststim.(['sub' subID])(trange_ttv,:);
            
            for c = 1:400
%                 currRT = allself.reac_time(startpoint+c-1);
%                 if currRT == -1
%                     currRT = 500;
%                 elseif currRT > 1599 
%                     currRT = 1599;
%                 end
                %if currRT > 400
%                     trange_ttv = ceil((currRT:(currRT+400))/2);
%                 else
%                     trange_ttv = 1:(currRT/2);
%                 end
                tmp2 = SlidingWindow(@nan_std,tmp(:,c),50,0.9);
                indvar(c) = mean(tmp2);
                %indvar(c) = SimplePLE(tmp);
            end
            
            indvar = indvar';
            indvar = double(indvar);
        elseif strcmpi(indvar_name,'amp_split')
            tmp = hilbert(prestim2.(['sub' subID]));
            indvar = mean(abs(tmp(950:1000,flatindices)),1);
            disp('Warning: hard-coded time')
            indvar = indvar';
            indvar = double(indvar);
        elseif strcmpi(indvar_name,'amp_split_win')
            tmp = hilbert(prestim2.(['sub' subID]));
            tmp = abs(tmp(trange_ttv,:));
            
            for c = 1:400
                tmp2 = SlidingWindow(@nan_std,tmp(:,c),50,0.9);
                indvar(c) = mean(tmp2);
                %indvar(c) = SimplePLE(tmp);
            end
            
            indvar = indvar';
            indvar = double(indvar);
        elseif strcmpi(indvar_name,'amp_slope')
            tmp = hilbert(poststim.(['sub' subID]));
            tmp = abs(tmp(trange_ttv,:));
            for c = 1:400
               tmp2 = polyfit(trange_ttv,tmp(:,c)',1);
               indvar(c) = tmp2(1);
            end
            indvar = indvar';
            indvar = double(indvar);
        end
        
        if EasyParse(varargin,'AccFilter','on')
%             rmtrials = [rmtrials find(isoutlier(indvar))'];
%             if ~strcmpi(depvar_name,'accuracy')
%                 rmtrials = [rmtrials find(isoutlier(depvar))'];
%             end
%             rmtrials = unique(rmtrials);
            
            indvar(rmtrials,:) = [];
            depvar(rmtrials,:) = [];
            rmtrials = [];
        end
        
        if ~EasyParse(varargin,'Subjectbased','false')
            indvar_record.(['sub' subID]) = indvar;
            depvar_record.(['sub' subID]) = depvar;
        else
            indvar_record{q} = indvar;
            depvar_record{q} = depvar;
        end

        
        if ~EasyParse(varargin,'Subjectbased','false')
            if strcmpi(depvar_name,'accuracy')
%                 if strcmpi(indvar_name,'phase')
%                     sinangle = sin(indvar);
%                     cosangle = cos(indvar);
%                     indvar = [sinangle cosangle];
%                 end
%                 [~,~,stats(q)] = mnrfit(indvar,ordinal(depvar));
                                     [r,p] = corr(indvar,depvar,'Type','Spearman');
                    stats(q).rho = r;
                    stats(q).p = p;
%                                              [r,p] = circ_corrcl(indvar,depvar);
%                                  stats(q).rho = r;
%                                  stats(q).p = p;
            elseif strcmpi(depvar_name,'acc_ttest')
                if ~CheckInput(varargin,'Bootstrap')
                x = indvar(find(depvar)); %accurate trials
                else
                    for c = 1:EasyParse(varargin,'Bootstrap')
                        
                    end
                end
                y = indvar(find(~depvar)); %inaccurate trials
                stats(q).p = ranksum(x,y,'tail','left');
                stats(q).meddif = median(y) - median(x);
                stats(q).meandif = mean(y) - mean(x);
            else
                if ~strcmpi(indvar_name,'phase') && ~strcmpi(indvar_name,'onset_phase') && ~strcmpi(depvar_name,'ersp_timecourse') && ~(iscell(electrode) && length(electrode) > 1)
                    [r,p] = corr(indvar,depvar,'Type','Spearman');
                    stats(q,1) = r;
                    stats(q,2) = p;
                elseif strcmpi(depvar_name,'ersp_timecourse') || (iscell(electrode) && length(electrode) > 1)
                    stats = NaN;
                else
                    [r,p] = circ_corrcl(indvar,depvar);
                    stats(q).rho = r;
                    stats(q).p = p;
                end
            end
            
        else
            %indvar = zscore(indvar);
            %indvar2 = vertcat(indvar2,indvar);
            %depvar2 = vertcat(depvar2,depvar);
        end
    end
    
    currentData = [];
end

indvar2 = zscore(indvar2);
if EasyParse(varargin,'Subjectbased','false') && ~strcmpi(depvar_name,'ersp_timecourse2')
    if strcmpi(depvar_name,'accuracy')
        rmindex = [];
       for c = 1:19
          if length(find(depvar_record{c} == 0)) < 5
              rmindex = [rmindex c];
          end
       end
       
       indvar_record(rmindex) = [];
       depvar_record(rmindex) = [];
    end
    tmp = cellfun(@transpose,indvar_record,'UniformOutput',false);
    tmp2 = double(cell2mat(tmp)');
    grouping = cellfun(@size,tmp,'UniformOutput',false); 
    grouping = cell2mat(cellfun(@transpose,grouping,'UniformOutput',false));
    grouping = grouping(2,:);
    tmp3 = [];
    for c = 1:length(depvar_record)
        tmp3 = cat(1,tmp3,c*ones(grouping(c),1));
    end
    
    tmp4 = double(cell2mat(cellfun(@transpose,depvar_record,'UniformOutput',false))');
    
    designTable = array2table([tmp2,tmp3,tmp4],'VariableNames',{indvar_name,'Participant',depvar_name});
    if ~strcmpi(depvar_name,'accuracy')
        stats = fitlme(designTable,[depvar_name ' ~ ' indvar_name ' + (1|Participant) + (-1 + ' indvar_name '|Participant)']);
    else
        stats = fitglme(designTable,[depvar_name ' ~ ' indvar_name ' + (1|Participant) + (-1 + ' indvar_name '|Participant)'],'Distribution','binomial','Link','logit');
    end
    
%     if strcmpi(depvar_name,'accuracy')
%         if strcmpi(indvar_name,'phase')
%             sinangle = sin(indvar2);
%             cosangle = cos(indvar2);
%             indvar2 = [sinangle cosangle];
%         end
%         [B,~,stats] = mnrfit(indvar2,ordinal(depvar2));
%         %                             [r,p] = circ_corrcl(indvar,depvar);
%         %                 stats(q).rho = r;
%         %                 stats(q).p = p;
%     elseif strcmpi(depvar_name,'RT')
%         if ~strcmpi(indvar_name,'phase')
%             [r,p] = corr(indvar2,depvar2,'Type','Spearman');
%             stats.rho = r;
%             stats.p = p;
%         else
%             [r,p] = circ_corrcl(indvar2,depvar2);
%             stats.rho = r;
%             stats.p = p;
%         end
%     end
elseif strcmpi(depvar_name,'ersp_timecourse2')
    stats = NaN;
end