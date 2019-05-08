function [ersp,pdifs,trange,pseudoersp] = NA_ERSP(trange,frange,conditions,electrode,splitmethod,varargin)

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

t0 = find(abs(timesTF) == min(abs(timesTF)));
if length(t0) > 1
   t0 = t0(1); 
end

for q = 1:length(fields)
    
    subID = extractAfter(fields{q},'sub');
    
    [~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(fields{q}));

    newindex = contains(filenames,extractAfter(fields{q},'sub'));
    
    newindex = fileindex(find(newindex));
    
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
                currentData = cat(4,currentData,tfdata(frange,:,:));
            catch
                disp(lasterror)
            end
        end
    end
    
    for i = 1:size(currentData,4)
        %first split the pseudotrials to get the different baselines
        if strcmpi(splitmethod,'power')
            pseudosplit = mean(mean(10*log10(abs(currentData(:,1:9,:,i).^2)),2),1);
            pseudosplit = squeeze(pseudosplit);
            
            ersp_split = pseudosplit > median(pseudosplit(flatindices{1}));
            ersp_split = cellstr(num2str(ersp_split));
            tmpflags = strcat(trialFlags.(fields{q}),ersp_split');
        elseif strcmpi(splitmethod,'cosphase')
            pseudosplit = cos(circ_mean(angle(currentData(:,1,:,i)),[],1));
            pseudosplit = squeeze(pseudosplit);
            
            ersp_split = pseudosplit > 0;
            ersp_split = cellstr(num2str(ersp_split));
            tmpflags = strcat(trialFlags.(fields{q}),ersp_split');
            
        elseif strcmpi(splitmethod,'concphase')
            allangles = wrapTo2Pi(circ_mean(angle(currentData(:,1,:,i)),[],1));
            
            cmean = circ_mean(allangles(flatindices{1}));
            cmean = wrapTo2Pi(cmean);
            
            concrange = [cmean-pi/2 cmean+pi/2];
            concrange = wrapTo2Pi(concrange);
            
            antirange = [cmean+pi/2 cmean+3*pi/2];
            antirange = wrapTo2Pi(antirange);
            
            tmp = zeros(1,400);
            if concrange(2) > concrange(1)
                tmp(intersect(find(allangles > concrange(1)),find(allangles < concrange(2)))) = 0;
            else
                tmp(intersect(find(allangles > concrange(1)),find(allangles < (2*pi)))) = 0;
                tmp(intersect(find(allangles < concrange(2)),find(allangles > 0))) = 0;
            end
            
            if antirange(2) > antirange(1)
                tmp(intersect(find(allangles > antirange(1)),find(allangles < antirange(2)))) = 1;
            else
                tmp(intersect(find(allangles > antirange(1)),find(allangles < (2*pi)))) = 1;
                tmp(intersect(find(allangles < antirange(2)),find(allangles > 0))) = 1;
            end
            
            mask = zeros(1,400);
            mask(flatindices{1}) = 1;
            tmp = tmp.*mask;
            
            tmp = tmp';
            tmp = cellstr(num2str(tmp));
            tmp = tmp';
            
            tmpflags = strcat(trialFlags.(fields{q}),tmp);
            
        end
        
        if ~strcmpi(splitmethod,'none')
            newconditions{1} = strcat(conditions{1},'0');
            newconditions{2} = strcat(conditions{2},'1');
            [flatindices2] = Ivar_findConds_3(newconditions,tmpflags);
        else
            flatindices2 = flatindices;
        end
        
        for c = 1:2
            pseudoersp{c}(:,:,q,i) = mean(10*log10(abs(currentData(:,10:(9+length(trange)),flatindices2{c},i).^2)),3);
        end
        
        %now split the real trials according to prestim power
        if strcmpi(splitmethod,'power')
            prestimsplit = mean(mean(10*log10(abs(currentData(:,(t0-8):t0,:,i).^2)),2),1);
            
            prestimsplit = squeeze(prestimsplit);
            
            ersp_split = prestimsplit > median(prestimsplit(flatindices{1}));
            ersp_split = cellstr(num2str(ersp_split));
            tmpflags = strcat(trialFlags.(fields{q}),ersp_split');
            
            [flatindices2] = Ivar_findConds_3(newconditions,tmpflags);
        elseif strcmpi(splitmethod,'cosphase')
            prestimsplit = cos(circ_mean(angle(currentData(:,t0,:,i).^2),[],1));
            
            prestimsplit = squeeze(prestimsplit);
            
            ersp_split = prestimsplit > 0;
            ersp_split = cellstr(num2str(ersp_split));
            tmpflags = strcat(trialFlags.(fields{q}),ersp_split');
            
            [flatindices2] = Ivar_findConds_3(newconditions,tmpflags);
        elseif strcmpi(splitmethod,'concphase')
            allangles = wrapTo2Pi(circ_mean(angle(currentData(:,t0,:,i)),[],1));
            
            cmean = circ_mean(allangles(flatindices{1}));
            cmean = wrapTo2Pi(cmean);
            
            concrange = [cmean-pi/2 cmean+pi/2];
            concrange = wrapTo2Pi(concrange);
            
            antirange = [cmean+pi/2 cmean+3*pi/2];
            antirange = wrapTo2Pi(antirange);
            
            tmp = zeros(1,400);
            if concrange(2) > concrange(1)
                tmp(intersect(find(allangles > concrange(1)),find(allangles < concrange(2)))) = 0;
            else
                tmp(intersect(find(allangles > concrange(1)),find(allangles < (2*pi)))) = 0;
                tmp(intersect(find(allangles < concrange(2)),find(allangles > 0))) = 0;
            end
            
            if antirange(2) > antirange(1)
                tmp(intersect(find(allangles > antirange(1)),find(allangles < antirange(2)))) = 1;
            else
                tmp(intersect(find(allangles > antirange(1)),find(allangles < (2*pi)))) = 1;
                tmp(intersect(find(allangles < antirange(2)),find(allangles > 0))) = 1;
            end
            
            mask = zeros(1,400);
            mask(flatindices{1}) = 1;
            tmp = tmp.*mask;
            
            tmp = tmp';
            tmp = cellstr(num2str(tmp));
            tmp = tmp';
            
            tmpflags = strcat(trialFlags.(fields{q}),tmp);
            
            [flatindices2] = Ivar_findConds_3(newconditions,tmpflags);
        elseif strcmpi(splitmethod,'none')
            flatindices2 = flatindices;
        end
        
        
        for c = 1:2
            ersp{c}(:,:,q,i) = mean(bsxfun(@minus,10*log10(abs(currentData(:,trange,flatindices2{c},i).^2)),mean(pseudoersp{c}(:,:,q,i),2)),3);
        end
        
        %normalize pseudotrial as well to compare
        for c = 1:2
            pseudoersp{c}(:,:,q,i) = bsxfun(@minus,pseudoersp{c}(:,:,q,i),mean(pseudoersp{c}(:,:,q,i),2));
        end
        
    end
    
    currentData = [];
end

for c = 1:2
    ersp{c} = squeeze(mean(ersp{c},4));
    pseudoersp{c} = squeeze(mean(pseudoersp{c},4));
end

for c = 1:length(trange)
    %pdifs(c) =  signrank(squeeze(mean(ersp{1}(:,c,:),1)),squeeze(mean(ersp{2}(:,c,:),1))));
    pdifs(c) = signrank(squeeze(mean(ersp{1}(:,c,:),1))-squeeze(mean(pseudoersp{1}(:,c,:),1)),squeeze(mean(ersp{2}(:,c,:),1)) - squeeze(mean(pseudoersp{2}(:,c,:),1)));
end

pdifs = mafdr(pdifs,'BHFDR',true);

if ~EasyParse(varargin,'Plot','off')
    NA_ERSPplot(ersp,orig_trange,trange,pdifs,'PlotPseudo',pseudoersp)
end

