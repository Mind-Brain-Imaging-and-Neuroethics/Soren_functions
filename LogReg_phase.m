function [stats] = LogReg_Phase(trange,frange,electrode,indvar,depvar)
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

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/allself.mat')


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
        
        
        accuracy = [];
        for c = 1:length(trialFlags.(['sub' subID]))
            if trialFlags.(['sub' subID]){c}(4) == '3'
                accuracy(c) = 1;
            elseif
                accuracy(c) = 0;
            end
        end
        
        
        
        currentData = mean(mean(currentData,1),2);
        currentData = squeeze(currentData);
        currentData = reshape(currentData,400,length(electrode));
        currentData = mean(currentData,2);
        
        %ersp = 10*log10(abs(tfdata).^2)
        
        %[rcorr(q) pcorr(q)] = circ_corrcl(angle(currentData),accuracy);
        sinangle = sin(angle(currentData));
        cosangle = cos(angle(currentData));
        
        [~,~,stats(q)] = mnrfit([sinangle cosangle],ordinal(accuracy),'Model','Ordinal');
    end
    
    
    currentData = [];
end

