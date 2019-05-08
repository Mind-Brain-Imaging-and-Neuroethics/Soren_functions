function [stats] = TrialBased_TF(trange,frange,electrode,indvar_name,depvar_name,varargin)
%logistic regression of prestimulus phase on accuracy

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/trialFlags_ivar.mat')
trialFlags = rmfield(trialFlags,'sub214');
trialFlags = rmfield(trialFlags,'sub224');

directory_name = '/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/raw_eeg/self/ICA_componentsrejected/All_eegcaps/';

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

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','TTV');

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TFparams.mat');

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

for q = 1:length(fields)
    
    newindex = contains(filenames,extractAfter(fields{q},'sub'));
    
    newindex = fileindex(find(newindex));
    
    if ~isempty(newindex)
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
                    currentData = cat(3,currentData,tfdata(frange,trange,:));
                catch
                    disp(lasterror)
                end
            end
        end
        
        currentData = mean(mean(currentData,1),2);
        currentData = squeeze(currentData);
        currentData = reshape(currentData,400,length(electrode));
        currentData = mean(currentData,2);
        
        rmtrials = [];
        depvar = [];
        indvar = [];
        
        if strcmpi(depvar_name,'accuracy')
            for c = 1:length(trialFlags.(['sub' subID]))
                if trialFlags.(['sub' subID]){c}(4) == '3'
                    depvar(c) = 1;
                else
                    depvar(c) = 0;
                end
%                 elseif trialFlags.(['sub' subID]){c}(4) == '2'
%                     depvar(c) = 0;
%                 else
%                     depvar(c) = NaN;
%                     rmtrials = [rmtrials c];
%                 end
            end
            depvar = depvar';
        elseif strcmpi(depvar_name,'RT')
            startPoint = find(allself.subid == str2num(subID),1);
            depvar = allself.reac_time(startPoint:startPoint+399);
            rmtrials = find(depvar == -1);
        end

        
        if strcmpi(indvar_name,'phase')
            %sinangle = sin(angle(currentData));
            %cosangle = cos(angle(currentData));
            indvar = angle(currentData);
        elseif strcmpi(indvar_name,'ersp')
            indvar = 10*log10(abs(currentData).^2); %not baseline corrected
%         elseif strcmpi(indvar_name,'ttv')
%             if orig_trange(1) < 0
%                indvar = mean(TTV.prestim.(['sub' subID])(:,1000+(orig_trange(1)/2):(1000+origtrange(2)/2)),1);
%                indvar = squeeze(indvar);
%                indvar = indvar';
%             end
        end
        
        indvar(rmtrials,:) = [];
        depvar(rmtrials,:) = [];
        
        if ~EasyParse(varargin,'Subjectbased','false')
            if strcmpi(depvar_name,'accuracy')
                if strcmpi(indvar_name,'phase')
                    sinangle = sin(indvar);
                    cosangle = cos(indvar);
                    indvar = [sinangle cosangle];
                end
                [~,~,stats(q)] = mnrfit(indvar,ordinal(depvar),'Model','Ordinal');
                %                             [r,p] = circ_corrcl(indvar,depvar);
                %                 stats(q).rho = r;
                %                 stats(q).p = p;
            elseif strcmpi(depvar_name,'RT')
                if ~strcmpi(indvar_name,'phase')
                    [r,p] = corr(indvar,depvar,'Type','Spearman');
                    stats(q).rho = r;
                    stats(q).p = p;
                else
                    [r,p] = circ_corrcl(indvar,depvar);
                    stats(q).rho = r;
                    stats(q).p = p;
                end
            end
            
        else
            indvar2 = vertcat(indvar2,indvar);
            depvar2 = vertcat(depvar2,depvar);
        end
    end
    
    
    currentData = [];
end

if EasyParse(varargin,'Subjectbased','false')
    if strcmpi(depvar_name,'accuracy')
        if strcmpi(indvar_name,'phase')
            sinangle = sin(indvar2);
            cosangle = cos(indvar2);
            indvar2 = [sinangle cosangle];
        end
        [~,~,stats] = mnrfit(indvar2,ordinal(depvar2),'Model','Ordinal');
        %                             [r,p] = circ_corrcl(indvar,depvar);
        %                 stats(q).rho = r;
        %                 stats(q).p = p;
    elseif strcmpi(depvar_name,'RT')
        if ~strcmpi(indvar_name,'phase')
            [r,p] = corr(indvar2,depvar2,'Type','Spearman');
            stats.rho = r;
            stats.p = p;
        else
            [r,p] = circ_corrcl(indvar2,depvar2);
            stats.rho = r;
            stats.p = p;
        end
    end
end