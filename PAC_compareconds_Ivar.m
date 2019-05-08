function [pac,pDifs,raw_pDifs] = PAC_compareconds_Ivar(conditions,electrode,Trange,Frange,varargin)

%loop for dividing by condition

load('/Volumes/Soren_SSD/Stuff_Ivar_SSD/eeg/allself.mat')
load('/Volumes/Soren_SSD/Stuff_Ivar_SSD/eeg/meandata_2.mat')
load('/Volumes/Soren_SSD/Stuff_Ivar_SSD/eeg/trialFlags_ivar_2.mat')
trialFlags = rmfield(trialFlags,'sub214');
trialFlags = rmfield(trialFlags,'sub224');
eeglab
directory_name = '/Volumes/Soren_SSD/Stuff_Ivar_SSD/eeg/raw_eeg/self/preprocessed_onlyContinuous/SelfICAMARA/';
cd(directory_name);

%%%%Let's select the relevant files, I am assuming here continuous resting-state files already clean
files = dir(fullfile(directory_name, '*.set'));

fileindex = find(~[files.isdir]);
for i = 1:length(fileindex)
    disp(['Now processing subject ' num2str(i)])
    
	filename = files(fileindex(i)).name;
	[PATH, NAME, EXT] = fileparts(filename);

	NAME = [NAME, EXT];
    
	EEG = pop_loadset( 'filename', filename, 'filepath', directory_name);
     
    subID = extractBetween(filename,'eeg_','_self');
    subID = subID{1};
    
    %subID = str2num(subID)
    
    %     startpos = find(allself.subid == subID,1)-1;
    %
    
    fnirscap = meandata.fnirscap(find(meandata.subid == str2num(subID)));
    electrodes = zeros(1,length(electrode));
    
    elecindex.fnirs = [13 14 19 20 40 51 45 46];
    
    elecindex.normal = [37 38 9 10 52 24 56 57];
    
    elecnames.fnirs = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};
    
    elecnames.normal = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};
    
    if fnirscap
        elec_num = elecindex.fnirs(find(contains(elecnames.fnirs,electrode)));
    else
        elec_num = elecindex.normal(find(contains(elecnames.normal,electrode)));
    end
    
    
    indx = [];
    for c = 1:length(EEG.event)
        if strcmpi(EEG.event(c).type,'S  5')
            indx = [indx c];
        end
    end
    
    for c = 1:400
        EEG.event(indx(c)).type = trialFlags.(['sub' subID]){c};
    end
    
    [~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(['sub' subID]));
    latencies = extractfield(EEG.event,'latency');
    
    indices{1} = latencies(indx(flatindices{1}));
    indices{2} = latencies(indx(flatindices{2}));
%     for q = 1:length(conditions)
%         for c = 1:400
%             if strcmpi(EEG.event(indx(c)).type,))
%                 indices{q} = [indices{q} EEG.event(c).latency];
%             end
%         end
%     end
    if ~CheckInput(varargin,'NormTrials')
        for c = 1:length(elec_num)
            [pac{1}(c,i,:),pac{2}(c,i,:),~,raw_pDifs(i,:)] = erpac_corr2(EEG.data(elec_num(c),:),EEG.srate,indices{1},indices{2},Trange,Frange(1),Frange(2),Frange(3),Frange(4));
        end
    else
%         for c = 1:EasyParse(varargin,'NormTrials')
%             [pac1(i,:),pac2(i,:),~,raw_pDifs(i,:)] = erpac_corr2(EEG.data(elec_num,:),EEG.srate,indices{1},indices{2},Trange,Frange(1),Frange(2),Frange(3),Frange(4));
%         end
        numTrials = cellfun('length',flatindices);
        minIndex = find(numTrials == min(numTrials));
        minTrials = min(numTrials);
        for c = 1:EasyParse(varargin,'NormTrials')
            bootIndices = randsample(indices{mod(minIndex,2)+1},minTrials);
            [pacBoot(i,c,:)] = erpac_corr(EEG.data(elec_num,:),EEG.srate,bootIndices,Trange,Frange(1),Frange(2),Frange(3),Frange(4));
        end
        [pac{minIndex}(i,:)] = erpac_corr(EEG.data(elec_num,:),EEG.srate,indices{minIndex},Trange,Frange(1),Frange(2),Frange(3),Frange(4));
    end
    %SelectTrials(EEG,subID,allself);
end

if length(electrode) > 1
    pac = cellfun(@mean,pac,'UniformOutput',false);
end
pac = cellfun(@squeeze,pac,'UniformOutput',false);

if ~CheckInput(varargin,'NormTrials')
    for c = 1:length(pac{1})
        pDifs(c) = signrank(pac{1}(:,c),pac{2}(:,c));
    end
else
    for cc = 1:length(pac{1})
        for c = 1:EasyParse(varargin,'NormTrials')
            tmpDifs(c,cc)= signrank(pac{minIndex}(:,c),pacBoot(:,c,cc));
        end
        %pDifs(cc) = EmpiricalBrownsMethod(horzcat(pacBoot(:,:,cc),pac{minIndex}(:,c)),tmpDifs(c,:));
        pDifs(cc) = Fisher_combine(tmpDifs(:,cc));
    end
end

if ~EasyParse(varargin,'MCorrect','off')
    if any(pDifs > 0.05) && any(pDifs < 0.05)
        pDifs = mafdr(pDifs);
    else
        pDifs = mafdr(pDifs,'BHFDR',true);
    end
end