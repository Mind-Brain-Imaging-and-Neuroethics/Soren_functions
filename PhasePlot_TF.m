function [testdifs,testsig,meanangles] = PhasePlot_TF(conditions,frange,trange,electrode,varargin)

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/trialFlags_ivar_2.mat')
trialFlags = rmfield(trialFlags,'sub214');
trialFlags = rmfield(trialFlags,'sub224');


directory_name = '/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/raw_eeg/self/ICA_componentsrejected/All_eegcaps/';

files = dir(fullfile(directory_name, '*.mat'));

fileindex = find(~[files.isdir]);

filenames = extractfield(files,'name');

currentData = [];

fields = fieldnames(trialFlags);

currentData = cell(1,length(conditions));

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

if ~EasyParse(varargin,'Subjectwise','true')
    for i = 1:length(fileindex)
        
        filename = files(fileindex(i)).name;
        [PATH, NAME, EXT] = fileparts(filename);
        
        NAME = [NAME, EXT];
        
        load([directory_name '/' filename])
        
        subID = extractBetween(filename,'tfdata_','_electrode');
        subID = subID{1};
        
        currElectrode = extractBetween(filename,'electrode','.mat');
        currElectrode = currElectrode{1};
        
        
        if CheckInput(varargin,'Electrode')
            electrode = EasyParse(varargin,'Electrode');
        else
            electrode = 'NA';
        end
        
        if any(contains(elecnames{str2num(currElectrode)},electrode)) || ~CheckInput(varargin,'Electrode')
            
            load([directory_name '/' filename])
            
            try
                [~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(['sub' subID]));
                for c = 1:length(flatindices)
                    currentData{c} = cat(3,currentData{c},tfdata(frange,trange,flatindices{c}));
                end
            catch
                disp(lasterror)
            end
        end
    end
    
    for c = 1:length(currentData)
        currentData{c} = mean(mean(currentData{c},1),2);
        currentData{c} = squeeze(currentData{c});
        currentData{c} = reshape(currentData{c},[],length(electrode));
        currentData{c} = mean(currentData{c},2);
    end
    
    figure
    for c = 1:length(flatindices)
        subplot(1,2,c)
        
        phaseplot = circ_plot(angle(reshape(currentData{c},[],1)),'hist');
        
        testsig(c) = circ_otest(angle(reshape(currentData{c},[],1)));
        meanangles(c) = circ_mean(angle(reshape(currentData{c},[],1)));
    end
    
    testdifs = circ_kuipertest(angle(reshape(currentData{1},[],1)),angle(reshape(currentData{2},[],1)));
    
else
    
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
                
                if CheckInput(varargin,'Electrode')
                    electrode = EasyParse(varargin,'Electrode');
                else
                    %electrode = horzcat(ele;
                end
                
                if any(contains(electrode,elecnames{str2num(currElectrode)}))
                    
                    load([directory_name '/' filename])
                    
                    try
                        [~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(['sub' subID]));
                        for c = 1:length(flatindices)
                            currentData{c} = cat(3,currentData{c},tfdata(frange,trange,flatindices{c}));
                        end
                    catch
                        disp(lasterror)
                    end
                end
            end
            
%             for c = 1:length(currentData)
%                currentData{c} = mean(mean(currentData{c},1),2);
%                currentData{c} = squeeze(currentData{c});
%                %currentData{c} = Mean_resample(currentData{c},length(electrode));
%                currentData{c} = reshape(currentData{c},length(flatindices{c}),length(electrode));
%                currentData{c} = mean(currentData{c},2);
%             end
            
            
            
            for c = 1:length(flatindices)
                try
                if c == 1
                   plotindex = (ceil(q/4)-1)*8+(q-(4*(ceil(q/4)-1)));
                elseif c == 2
                    plotindex = (ceil(q/4)-1)*8+(q-(4*(ceil(q/4)-1)))+4;
                end
                %plotindex = ((q-1)*2+c);
                subplot(5,8,plotindex)
                angles{c} = angle(currentData{c});
                angles{c} = squeeze(circ_mean(circ_mean(angles{c},[],2),[],1));
                
                phaseplot = circ_plot(angles{c},'hist');
                
                testsig(c,q) = circ_otest(angles{c});
                meanangles(c,q) = circ_mean(angles{c});
                
                catch
                    disp(['Subject ' subID])
                    disp(lasterror)
                end
            end
            %testdifs(q) = circ_kuipertest(angles{1},angles{2});

            
            
        end
        
        currentData = cell(1,length(conditions));
    end
end