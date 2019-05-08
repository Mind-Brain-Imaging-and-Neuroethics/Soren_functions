function [stats,indvar_subs,depvar_subs] = SubjectCorr_Ivar(trange,frange,electrode,indvar_name,depvar_name,varargin)
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

%load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','TTV');

if ~EasyParse(varargin,'LongEpochs','true')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TFparams.mat');
else
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TFparams_long.mat');
end

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/corrRest12.mat','allMeasures','EEGM');

if strcmpi(indvar_name,'ple_elec')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/corrRest12.mat','PLE_JF_onlysubs');
    rest_elec = PLE_JF_onlysubs;
elseif strcmpi(indvar_name,'dfa_elec')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/corrRest12.mat','DFA_onlysubs');
    rest_elec = DFA_onlysubs;
elseif strcmpi(indvar_name,'ple_irasa_elec')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/corrRest12.mat','PLE_IRASA_onlysubs');
    rest_elec = PLE_IRASA_onlysubs;
elseif strcmpi(indvar_name,'oscifrac_elec')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/corrRest12.mat','OsciFrac_onlysubs');
    rest_elec = OsciFrac_onlysubs;
elseif strcmpi(indvar_name,'dfa_alpha_elec')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/corrRest12.mat','DFA_amp_clean_alpha');
    DFA_amp_clean_alpha = GetSubsReal(DFA_amp_clean_alpha,EEGM,meandata);
    rest_elec = DFA_amp_clean_alpha;
end

if strcmpi(indvar_name,'bandpower') || strcmpi(depvar_name,'bandpower')
    bpargs = EasyParse(varargin,'BPargs');
    if strcmpi(bpargs{1},'rel')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/Calculated_measures/all_rel_bandpower.mat','outmeasures')
    elseif strcmpi(bpargs{1},'abs')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/Calculated_measures/all_bandpower.mat','outmeasures')
    elseif strcmpi(bpargs{1},'rel_osci')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/Calculated_measures/Rel_osci_bandpower.mat','outmeasures')
    elseif strcmpi(bpargs{1},'abs_osci')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/Calculated_measures/Osci_bandpower_ple.mat','outmeasures')
    end
    bandwidths = {'delta','theta','alpha','beta','gamma'};
    rest_elec2 = outmeasures(:,:,find(contains(bandwidths,bpargs{2})));
    rest_elec2 = GetSubsReal_34(rest_elec2,EEGM,meandata);
end

if strcmpi(indvar_name,'bpcustom') || strcmpi(depvar_name,'bpcustom')
    if EasyParse(varargin,'BPmode','rel')
        outmeasures = Loop_measure({@(EEG)Rel_Bandpower_EEG_wrapper(EEG,frange)},[],'Dir','/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/raw_eeg/all_rest/2_7minclean/RestICAMARA/');
    else
        outmeasures = Loop_measure({@(EEG)Bandpower_EEG_wrapper(EEG,frange)},[],'Dir','/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/raw_eeg/all_rest/2_7minclean/RestICAMARA/');
    end
    
    rest_elec2 = GetSubsReal_34(outmeasures,EEGM,meandata);
end

if strcmpi(indvar_name,'ple_hf') || strcmpi(depvar_name,'ple_hf') || strcmpi(indvar_name,'ple_lf') || strcmpi(depvar_name,'ple_lf')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/Calculated_measures/Osci_bandpower_ple.mat','outmeasures')
    if strcmpi(indvar_name,'ple_hf') || strcmpi(depvar_name,'ple_hf')
        rest_elec3 = GetSubsReal_34(outmeasures(:,:,7),EEGM,meandata);
    else
        rest_elec3 = GetSubsReal_34(outmeasures(:,:,6),EEGM,meandata);
    end
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

if strcmpi(indvar_name,'ttv') || strcmpi(depvar_name,'ttv')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','poststim')
    %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat','prestim2')
    %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','prestim')
    %poststim = prestim2;
    %clear prestim
    if CheckInput(varargin,'Conditions')
        [~,datattv,tpointsTTV] = TTV_compareconds(poststim,electrode,EasyParse(varargin,'Conditions'),'Plot','off','MCorrect','off','SigDifs','off');
    else
        [~,datattv,tpointsTTV] = TTV_compareconds(poststim,electrode,{'xxxxx','43210'},'Plot','off','MCorrect','off','SigDifs','off');
    end
    
    [~,t1] = min(abs(tpointsTTV-orig_trange(1)));
    [~,t2] = min(abs(tpointsTTV-orig_trange(2)));
    ttv_trange = t1:t2;
    disp(['Taking time points ' num2str(t1) ' to ' num2str(t2) '...'])
    datattv = datattv{1};
    %clear poststim;
end
if strcmpi(indvar_name,'ttv_window') || strcmpi(depvar_name,'ttv_window')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','poststim')
    if CheckInput(varargin,'Conditions')
        [~,datattv,tpointsTTV] = TTV_compareconds(poststim,electrode,EasyParse(varargin,'Conditions'),'Plot','off','MCorrect','off','SigDifs','off','SlidingWindow',[50 0.9]);
    else
        [~,datattv,tpointsTTV] = TTV_compareconds(poststim,electrode,{'xxxxx','43210'},'Plot','off','MCorrect','off','SigDifs','off','SlidingWindow',[50 0.9]);
    end
    [~,t1] = min(abs(tpointsTTV-orig_trange(1)));
    [~,t2] = min(abs(tpointsTTV-orig_trange(2)));
    ttv_trange = t1:t2;
    disp(['Taking time points ' num2str(t1) ' to ' num2str(t2) '...'])
    datattv = datattv{1};
    %clear poststim;
end
if strcmpi(indvar_name,'nonadd_index') || strcmpi(depvar_name,'nonadd_index')
    if ~CheckInput(varargin,'TTVband')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','poststim')
    else
        load(['/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc_' EasyParse(varargin,'TTVband') '.mat'],['poststim_' EasyParse(varargin,'TTVband')]);
        if EasyParse(varargin,'TTVband','alpha')
            poststim = poststim_alpha;
            clear poststim_alpha;
        elseif EasyParse(varargin,'TTVband','beta')
            poststim = poststim_beta;
            clear poststim_beta;
        end
    end
    %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat','prestim2')
    if CheckInput(varargin,'Conditions')
        [~,datattv] = TTV_compareconds(poststim,electrode,EasyParse(varargin,'Conditions'),'Plot','off','MCorrect','off','SigDifs','off','PrestimSplit',{'Trange',[1900 2000],'Frange',[0 3],'Criterion','amp_split'});
    else
        [~,datattv] = TTV_compareconds(poststim,electrode,{'1xxxx','1xxxx'},'Plot','off','MCorrect','off','SigDifs','off','PrestimSplit',{'Trange',[1900 2000],'Frange',[0 3],'Criterion','amp_split'});
    end
    datattv1 = datattv{1};
    datattv2 = datattv{2};
end

if strcmpi(indvar_name,'nonaddindex_zirui') || strcmpi(depvar_name,'nonaddindex_zirui')
      if ~CheckInput(varargin,'TTVband')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','poststim')
    else
        load(['/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc_' EasyParse(varargin,'TTVband') '.mat'],['poststim_' EasyParse(varargin,'TTVband')]);
        if EasyParse(varargin,'TTVband','alpha')
            poststim = poststim_alpha;
            clear poststim_alpha;
        elseif EasyParse(varargin,'TTVband','beta')
            poststim = poststim_beta;
            clear poststim_beta;
        end
    end
    %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat','prestim2')
    if CheckInput(varargin,'Conditions')
        %[~,datattvna] = TTV_compareconds(poststim,electrode,EasyParse(varargin,'Conditions'),'Plot','off','MCorrect','off','SigDifs','off','NormRange',[1 1]);
        [~,datattvna] = TTV_compareconds(poststim,electrode,EasyParse(varargin,'Conditions'),'Plot','off','MCorrect','off','SigDifs','off','Pseudotrial','broadband','Limits',[1 500]);
    else
        [~,datattvna] = TTV_compareconds(poststim,electrode,{'1xxxx','43210'},'Plot','off','MCorrect','off','SigDifs','off','Pseudotrial','broadband','Limits',[1 500]);
    end
    warning('Pseudotrial normalization is on')
    %[~,datattv2] = TTV_compareconds(prestim2,electrode,{'1xxxx','43210'},'Plot','off','MCorrect','off','SigDifs','off');
    datattvna = datattvna{1};
    %datattv2 = datattv2{1};
end

if strcmpi(indvar_name,'pseudo_nonadd') || strcmpi(depvar_name,'pseudo_nonadd')
    if ~CheckInput(varargin,'TTVband')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','poststim')
    else
        load(['/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc_' EasyParse(varargin,'TTVband') '.mat'],['poststim_' EasyParse(varargin,'TTVband')]);
        if EasyParse(varargin,'TTVband','alpha')
            poststim = poststim_alpha;
            clear poststim_alpha;
        elseif EasyParse(varargin,'TTVband','beta')
            poststim = poststim_beta;
            clear poststim_beta;
        end
    end
    %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat','prestim2')
    if CheckInput(varargin,'Conditions')
        [~,datattvp] = TTV_compareconds(poststim,electrode,EasyParse(varargin,'Conditions'),'Plot','off','MCorrect','off','SigDifs','off','PrestimSplit',{'Trange',[1900 2000],'Frange',[0 3],'Criterion','amp_split'},'PseudotrialSplit','broadband','Limits',[1 450]);
    else
        [~,datattvp] = TTV_compareconds(poststim,electrode,{'1xxxx','1xxxx'},'Plot','off','MCorrect','off','SigDifs','off','PrestimSplit',{'Trange',[1900 2000],'Frange',[0 3],'Criterion','amp_split'},'PseudotrialSplit','broadband','Limits',[1 450]);
    end
    %datattv1 = datattv{1};
    %datattv2 = datattv{2};
end

if strcmpi(indvar_name,'nonadd_slope') || strcmpi(depvar_name,'nonadd_slope')
    if ~CheckInput(varargin,'TTVband')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','poststim')
    else
        load(['/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc_' EasyParse(varargin,'TTVband') '.mat'],['poststim_' EasyParse(varargin,'TTVband')]);
        if EasyParse(varargin,'TTVband','alpha')
            poststim = poststim_alpha;
            clear poststim_alpha;
        elseif EasyParse(varargin,'TTVband','beta')
            poststim = poststim_beta;
            clear poststim_beta;
        end
    end
    %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat','prestim2')
    if CheckInput(varargin,'Conditions')
        [~,datattvsl] = TTV_compareconds(poststim,electrode,EasyParse(varargin,'Conditions'),'Plot','off','MCorrect','off','SigDifs','off','NormRange',[1 1]);
    else
        [~,datattvsl] = TTV_compareconds(poststim,electrode,{'1xxxx','43210'},'Plot','off','MCorrect','off','SigDifs','off');
    end
    warning('Trial normalization is on')
    %[~,datattv2] = TTV_compareconds(prestim2,electrode,{'1xxxx','43210'},'Plot','off','MCorrect','off','SigDifs','off');
    datattvsl = datattvsl{1};
    %datattv2 = datattv2{1};
end
if strcmpi(indvar_name,'interstim_slope') || strcmpi(depvar_name,'interstim_slope')
    if ~CheckInput(varargin,'TTVband')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','interstim')
    else
        load(['/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc_' EasyParse(varargin,'TTVband') '.mat'],['interstim_' EasyParse(varargin,'TTVband')]);
        if EasyParse(varargin,'TTVband','alpha')
            interstim = interstim_alpha;
            clear interstim_alpha;
        elseif EasyParse(varargin,'TTVband','beta')
            interstim = interstim_beta;
            clear interstim_beta;
        end
    end
    %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat','prestim2')
    if CheckInput(varargin,'Conditions')
        [~,datattvsl_inter] = TTV_compareconds(interstim,electrode,EasyParse(varargin,'Conditions'),'Plot','off','MCorrect','off','SigDifs','off');
    else
        [~,datattvsl_inter] = TTV_compareconds(interstim,electrode,{'4xxxx','43210'},'Plot','off','MCorrect','off','SigDifs','off');
    end
    %[~,datattv2] = TTV_compareconds(prestim2,electrode,{'1xxxx','43210'},'Plot','off','MCorrect','off','SigDifs','off');
    datattvsl_inter = datattvsl_inter{1};
    if CheckInput(varargin,'Conditions')
        condns = EasyParse(varargin,'Conditions');
        if strcmpi(condns{1},'3xxxx')
            datattvsl_inter = datattvsl_inter(:,1:125,:);
        elseif strcmpi(condns{1},'2xxxx')
            datattvsl_inter = datattvsl_inter(:,1:84,:);
        end
    end
    %datattv2 = datattv2{1};
end    
if strcmpi(indvar_name,'ttv_prestim')
    if ~CheckInput(varargin,'TTVband')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','prestim')
    else
        load(['/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc_' EasyParse(varargin,'TTVband') '.mat'],['prestim_' EasyParse(varargin,'TTVband')]);
        
        if EasyParse(varargin,'TTVband','alpha')
            prestim = prestim_alpha;
            clear prestim_alpha;
        elseif EasyParse(varargin,'TTVband','beta')
            prestim = prestim_beta;
            clear prestim_beta;
        end
    end
    %load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat','prestim2')
    if CheckInput(varargin,'Conditions')
        [~,datattv_pre,tpointsTTV] = TTV_compareconds(prestim,electrode,EasyParse(varargin,'Conditions'),'Plot','off','MCorrect','off','SigDifs','off');
    else
        [~,datattv_pre,tpointsTTV] = TTV_compareconds(prestim,electrode,{'xxxxx','43210'},'Plot','off','MCorrect','off','SigDifs','off');
    end
    [~,t1] = min(abs(tpointsTTV-orig_trange(1)));
    [~,t2] = min(abs(tpointsTTV-orig_trange(2)));
    ttv_trange = t1:t2;
    disp(['Taking time points ' num2str(t1) ' to ' num2str(t2) '...'])
    datattv_pre = datattv_pre{1};
end

if strcmpi(indvar_name,'p300') || strcmpi(depvar_name,'p300') || strcmpi(indvar_name,'n200') || strcmpi(depvar_name,'n200') || strcmpi(indvar_name,'p200') || strcmpi(depvar_name,'p200') || strcmpi(indvar_name,'n100') || strcmpi(depvar_name,'n100')
    if EasyParse(varargin,'ERPstim',1)
            load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','interstim')
        poststim = interstim;
        clear interstim;
    else
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','poststim')
    end
end

if strcmpi(indvar_name,'MI') || strcmpi(depvar_name,'MI') || strcmpi(indvar_name,'MIdifs') || strcmpi(depvar_name,'MIdifs')
   load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/MItestvars.mat','MIvals','pseudoMIvals')
end

if strcmpi(indvar_name,'erpac') || strcmpi(depvar_name,'erpac')
    erpac_frange = EasyParse(varargin,'Erpac_frange');
    erpac_trange = EasyParse(varargin,'Erpac_trange');
    if CheckInput(varargin,'Conditions')
        [pac] = PAC_compareconds_Ivar(EasyParse(varargin,'Conditions'),electrode,erpac_trange,erpac_frange,'MCorrect','off');
    else
        [pac] = PAC_compareconds_Ivar({'xxxxx','xxxxx'},electrode,erpac_trange,erpac_frange,'MCorrect','off');
    end
    pac = pac{1};
end

savedata = [];
for q = 1:length(fields)
    
    subID = extractAfter(fields{q},'sub');
    
    newindex = contains(filenames,subID);
    
    if CheckInput(varargin,'Conditions')
        [~,~,flatindices] = Ivar_findConds_2(EasyParse(varargin,'Conditions'),trialFlags.(fields{q}));
        flatindices = flatindices{1};
    else
        flatindices = 1:400;
    end
    
    newindex = fileindex(find(newindex));
    
%     if strcmpi(indvar_name,'ersp')
%         filename = files(fileindex(i)).name;
%         
%         EEG = pop_loadset('filename', filename, 'filepath', directory_name);
%         
%             
%         if strcmpi(EEG.chanlocs(1).labels,'Fp1')
%             EEG.data = EEG.data([37 38 9 10 52 24 56 57],:,:); %normal EEG cap
%             EEG.chanlocs = EEG.chanlocs([37 38 9 10 52 24 56 57]);
%         else
%             EEG.data = EEG.data([13 14 19 20 40 51 45 46],:,:); %fnirs EEG cap
%             EEG.chanlocs = EEG.chanlocs([13 14 19 20 40 51 45 46]);
%         end
%         EEG.nbchan = 8;
% 
%     end

    if strcmpi(indvar_name,'ersp') || strcmpi(depvar_name,'ersp') || strcmpi(indvar_name,'phase') || strcmpi(depvar_name,'phase') || strcmpi(indvar_name,'std_ersp') || strcmpi(indvar_name,'itc') || strcmpi(depvar_name,'itc') || strcmpi(depvar_name,'raw_ersp')
        if ~isempty(newindex)
            for i = 1:length(newindex)
                
                filename = files(newindex(i)).name;
                [PATH, NAME, EXT] = fileparts(filename);
                
                NAME = [NAME, EXT];
                
                %             subID = extractBetween(filename,'tfdata_','_electrode');
                %             subID = subID{1};
                
                currElectrode = extractBetween(filename,'electrode','.mat');
                currElectrode = currElectrode{1};
                
                if any(contains(electrode,elecnames{str2num(currElectrode)}))
                    
                    load([directory_name '/' filename])
                    
                    try
                        currentData = cat(3,currentData,tfdata(frange,:,flatindices));
                        warning('Multi-electrode ERSP is not fixed yet')
                    catch
                        disp(lasterror)
                    end
                end
            end
            
            
%             currentData = mean(mean(currentData,1),2);
%             currentData = squeeze(currentData);
%             currentData = reshape(currentData,400,length(electrode));
%             currentData = mean(currentData,2);
        end
    end
    rmtrials = [];
    indvar = [];
    depvar = [];
    
    if strcmpi(indvar_name,'p300') || strcmpi(depvar_name,'p300') || strcmpi(indvar_name,'n200') || strcmpi(depvar_name,'n200') || strcmpi(indvar_name,'p200') || strcmpi(depvar_name,'p200') || strcmpi(indvar_name,'n100') || strcmpi(depvar_name,'n100')
        if iscell(electrode)
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
        else
            poststim.(['sub' subID]) = squeeze(mean(poststim.(['sub' subID])(electrode,:,:),1));
        end
    end
    
    if iscell(electrode)
        fnirscap = meandata.fnirscap(find(meandata.subid == str2num(subID)));
        ttv_electrodes = zeros(1,length(electrode));
        
        ttv_elecindex.fnirs = [13 14 19 20 40 51 45 46];
        
        ttv_elecindex.normal = [37 38 9 10 52 24 56 57];
        
        ttv_elecnames.fnirs = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};
        
        ttv_elecnames.normal = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};
        
        if fnirscap
            subElecs = ttv_elecindex.fnirs(find(contains(ttv_elecnames.fnirs,electrode)));
        else
            subElecs = ttv_elecindex.normal(find(contains(ttv_elecnames.normal,electrode)));
        end
    else
        subElecs = electrode;
    end
    
    if strcmpi(depvar_name,'meanacc')
        %             for c = 1:length(trialFlags.(['sub' subID]))
        %                 if trialFlags.(['sub' subID]){c}(4) == '3'
        %                     depvar(c) = 1;
        %                 else
        %                     depvar(c) = 0;
        %                 end
        % %                 elseif trialFlags.(['sub' subID]){c}(4) == '2'
        % %                     depvar(c) = 0;
        % %                 else
        % %                     depvar(c) = NaN;
        % %                     rmtrials = [rmtrials c];
        % %                 end
        %             end
        %             depvar = depvar';
        depvar = meandata.MeanAcc(find(meandata.subid == str2num(subID)));
    elseif strcmpi(depvar_name,'accuracy')
        currBehavData = allself(find(allself.subid == str2num(subID)),:);
        currBehavData = currBehavData.accuracy;
        currBehavData(find(currBehavData == -1)) = 0;
        %[~,~,flatindices] = Ivar_findConds_2({'1xxxx','xxxxx'},trialFlags.(fields{q}));
        %flatindices = flatindices{1};
        depvar = mean(currBehavData(flatindices));
    elseif strcmpi(depvar_name,'RT')
                    startPoint = find(allself.subid == str2num(subID),1);
                    depvar = allself.reac_time(startPoint:startPoint+399);
                    depvar = depvar(flatindices);
                    depvar(find(depvar == -1)) = [];
        %depvar = meandata.MeanRT(find(meandata.subid == str2num(subID)));
    elseif strcmpi(depvar_name,'EfScore')
        depvar = meandata.MeanRT(find(meandata.subid == str2num(subID)))/meandata.MeanAcc(find(meandata.subid == str2num(subID)));
    elseif strcmpi(depvar_name,'ttv')
        depvar = mean(mean(datattv(:,trange,q)));
    elseif strcmpi(depvar_name,'ttv_window')
        depvar = mean(mean(datattv(:,ttv_trange,q)));
    elseif strcmpi(depvar_name,'phase')
        depvar = circ_mean(angle(currentData));
    elseif strcmpi(depvar_name,'ersp')
        depvar = mean(mean(mean(bsxfun(@minus,10*log10(abs(currentData(:,trange,:)).^2),mean(mean(10*log10(abs(currentData(:,1:60,:)).^2),3),2))))); 
    elseif strcmpi(depvar_name,'raw_ersp')
        depvar = mean(mean(mean(10*log10(abs(currentData(:,trange,:)).^2)))); 
    elseif strcmpi(depvar_name,'itc')
        depvar = mean(mean(abs(mean((currentData(:,trange,:)./abs(currentData(:,trange,:))),3))));
        depvar = abs(depvar);
    elseif strcmpi(depvar_name,'nonadd_slope')
        slope_post = polyfit(1:200,squeeze(mean(datattvsl(:,1:200,q),1)),1);
        slope_post = slope_post(1);
        %slope_pre = polyfit(1:200,squeeze(mean(datattv2(:,801:1000,q),1)),1);
        %slope_pre = slope_pre(1);
        %depvar = slope_pre-slope_post;
        depvar = slope_post;
    elseif strcmpi(depvar_name,'interstim_slope')
        slope_inter = polyfit(1:size(datattvsl_inter,2),squeeze(mean(datattvsl_inter(:,:,q),1)),1);
        depvar = slope_inter(1);
    elseif strcmpi(depvar_name,'nonadd_index')
        depvar = trapz(mean(datattv2(:,1:200,q),1))-trapz(mean(datattv1(:,1:200,q),1));
    elseif strcmpi(depvar_name,'nonaddindex_zirui')
        depvar = trapz(mean(datattvna(:,1:200,q),1));
    elseif strcmpi(depvar_name,'pseudo_nonadd')
        depvar = trapz(mean(datattvp{1}(:,1:200,q),1)) - trapz(mean(datattvp{2}(:,1:200,q),1));
    elseif strcmpi(depvar_name,'amp_win')
        tmp = hilbert(squeeze(mean(poststim.(['sub' subID]),1)));
        tmp = abs(tmp(ceil((orig_trange(1)/2)):ceil(orig_trange(2)/2),:));
        %tmp = poststim.(['sub' subID])(trange_ttv,:);
        
        for c = 1:400
            tmp2 = SlidingWindow(@nan_std,tmp(:,c),50,0.9);
            depvar(c) = mean(tmp2);
            %indvar(c) = SimplePLE(tmp);
        end
        
        depvar = depvar';
        depvar = double(depvar);
    elseif strcmpi(depvar_name,'p300')
        %[~,~,flatindices] = Ivar_findConds_2({'1xxxx','1xxxx'},trialFlags.(fields{q}));
        tmp = mean(poststim.(fields{q})(125:225,flatindices),2);
        pks = findpeaks(tmp);
        depvar = max(tmp(pks.loc));
    elseif strcmpi(depvar_name,'n200')
        [~,~,flatindices] = Ivar_findConds_2({'1xxxx','1xxxx'},trialFlags.(fields{q}));
        tmp = mean(poststim.(fields{q})(50:150,flatindices{1}),2);
        tmp = -tmp;
        pks = findpeaks(tmp);
        depvar = max(tmp(pks.loc));
        depvar = -depvar;
    elseif strcmpi(depvar_name,'p200')
        [~,~,flatindices] = Ivar_findConds_2({'1xxxx','1xxxx'},trialFlags.(fields{q}));
        tmp = mean(poststim.(fields{q})(75:150,flatindices{1}),2);
        pks = findpeaks(tmp);
        depvar = max(tmp(pks.loc));
    elseif strcmpi(depvar_name,'n100')
        [~,~,flatindices] = Ivar_findConds_2({'1xxxx','1xxxx'},trialFlags.(fields{q}));
        tmp = mean(poststim.(fields{q})(40:80,flatindices{1}),2);
        tmp = -tmp;
        pks = findpeaks(tmp);
        depvar = max(tmp(pks.loc));
        depvar = -depvar;
    elseif strcmpi(depvar_name,'self_RT_1')
        depvar = meandata.SelfMeanRT1(find(meandata.subid == str2num(subID)));
    elseif strcmpi(depvar_name,'MI')
        delayval = EasyParse(varargin,'Delay');
        depvar = mean(mean(MIvals{q,delayval}(:,subElecs)));
    elseif strcmpi(depvar_name,'MIdifs')
        delayval = EasyParse(varargin,'Delay');
        depvar = mean(mean(MIvals{q,delayval}(:,subElecs))) - mean(mean(pseudoMIvals{q,delayval}(:,subElecs)));
    elseif strcmpi(depvar_name,'erpac')
        depvar = mean(pac(q,:));
        %elseif strcmpi(depvar_name,'
    elseif strcmpi(depvar_name,'accslope')
        
    elseif strcmpi(depvar_name,'bandpower') || strcmpi(depvar_name,'bpcustom')
         depvar = squeeze(mean(rest_elec2(find(meandata.subid == str2num(subID)),subElecs)));
    elseif strcmpi(depvar_name,'ple_lf') || strcmpi(depvar_name,'ple_hf')
        depvar = squeeze(mean(rest_elec3(find(meandata.subid == str2num(subID)),subElecs)));
    end
    
    %bsxfun(@plus,10*log10(abs(tfdata(:,:,c)).^2)',-powbase)
    if strcmpi(indvar_name,'phase')
        indvar = circ_mean(angle(currentData));
    elseif strcmpi(indvar_name,'ersp')
        indvar = mean(mean(mean(bsxfun(@minus,10*log10(abs(currentData(:,trange,:)).^2),mean(mean(10*log10(abs(currentData(:,1:60,:)).^2),3),2))))); 
    elseif strcmpi(indvar_name,'std_ersp')
        indvar = 10*log10(abs(currentData).^2); %not baseline corrected
        indvar = std(indvar);
        %         elseif strcmpi(indvar_name,'ttv')
        %             if orig_trange(1) < 0
        %                indvar = mean(TTV.prestim.(['sub' subID])(:,1000+(orig_trange(1)/2):(1000+origtrange(2)/2)),1);
        %                indvar = squeeze(indvar);
        %                indvar = indvar';
        %             end
    elseif strcmpi(indvar_name,'itc')
        %savedata = cat(3,savedata,currentData);
        indvar = mean(mean(abs(mean((currentData(:,trange,:)./abs(currentData(:,trange,:))),3))));
        indvar = abs(indvar);
    elseif strcmpi(indvar_name,'itc_timecourse')
        indvar = mean(abs(mean((currentData(:,trange,:)./abs(currentData(:,trange,:))),3)),1);
        indvar = abs(indvar);
        indvar = squeeze(indvar);
    elseif strcmpi(indvar_name,'ttv')
        indvar = mean(mean(datattv(:,ttv_trange,q)));
    elseif strcmpi(indvar_name,'ttv_window')
        indvar = mean(mean(datattv(:,ttv_trange,q)));
    %elseif strcmpi(indvar_name,'amp_win')
        %indvar = 
    elseif strcmpi(indvar_name,'nonadd_slope')
        slope_post = polyfit(1:200,squeeze(mean(datattvsl(:,1:200,q),1)),1);
        slope_post = slope_post(1);
        %slope_pre = polyfit(1:200,squeeze(mean(datattv2(:,801:1000,q),1)),1);
        %slope_pre = slope_pre(1);
        %indvar = slope_pre-slope_post;
        indvar = slope_post;
    elseif strcmpi(indvar_name,'pseudo_nonadd')
        indvar = trapz(mean(datattvp{1}(:,1:200,q),1)) - trapz(mean(datattvp{2}(:,1:200,q),1));
    elseif strcmpi(indvar_name,'interstim_slope')
        slope_inter = polyfit(1:size(datattvsl_inter,2),squeeze(mean(datattvsl_inter(:,:,q),1)),1);
        indvar = slope_inter(1);
    elseif strcmpi(indvar_name,'ttv_prestim')
        indvar = mean(mean(datattv_pre(:,ttv_trange,q)));
    elseif strcmpi(indvar_name,'amplitude')
         tmp = hilbert(squeeze(mean(poststim.(['sub' subID]),1)));
        tmp = abs(tmp(ceil((orig_trange(1)/2)):ceil(orig_trange(2)/2),:));
        tmp = squeeze(mean(tmp,1));
        indvar = tmp;
        indvar = indvar';
        indvar = double(indvar);
    elseif strcmpi(indvar_name,'amp_win')
        tmp = hilbert(squeeze(mean(poststim.(['sub' subID]),1)));
        tmp = abs(tmp(ceil((orig_trange(1)/2)):ceil(orig_trange(2)/2),:));
        %tmp = poststim.(['sub' subID])(trange_ttv,:);
        
        for c = 1:400
            tmp2 = SlidingWindow(@nan_std,tmp(:,c),50,0.9);
            indvar(c) = mean(tmp2);
            %indvar(c) = SimplePLE(tmp);
        end
        
        indvar = indvar';
        indvar = double(indvar);
    elseif strcmpi(indvar_name,'nonadd_index')
        indvar = trapz(mean(datattv2(:,1:200,q),1))-trapz(mean(datattv1(:,1:200,q),1));
    elseif strcmpi(indvar_name,'ple')
        indvar = allMeasures.PLE_raw_clean(find(meandata.subid == str2num(subID)));
    elseif strcmpi(indvar_name,'dfa')
        indvar = allMeasures.DFA_amp_clean(find(meandata.subid == str2num(subID)));
    elseif strcmpi(indvar_name,'p300')
        [~,~,flatindices] = Ivar_findConds_2({'1xxxx','1xxxx'},trialFlags.(fields{q}));
        tmp = mean(poststim.(fields{q})(125:225,flatindices{1}),2);
        pks = findpeaks(tmp);
        indvar = max(tmp(pks.loc));
    elseif strcmpi(indvar_name,'erpac')
        indvar = mean(pac(q,:));
    elseif strcmpi(indvar_name,'bandpower') || strcmpi(indvar_name,'bpcustom')
        indvar = squeeze(mean(rest_elec2(find(meandata.subid == str2num(subID)),subElecs)));
    elseif strcmpi(indvar_name,'ple_elec') || strcmpi(indvar_name,'dfa_elec') || strcmpi(indvar_name,'ple_irasa_elec') || strcmpi(indvar_name,'oscifrac_elec') || strcmpi(indvar_name,'dfa_alpha_elec')
        indvar = squeeze(mean(rest_elec(find(meandata.subid == str2num(subID)),subElecs)));
    elseif strcmpi(indvar_name,'ple_lf') || strcmpi(indvar_name,'ple_hf')
        indvar = squeeze(mean(rest_elec3(find(meandata.subid == str2num(subID)),subElecs)));
    end
    
    %indvar_subs(rmtrials,:) = [];
    %depvar(rmtrials,:) = [];
    if ~strcmpi(indvar_name,'itc_timecourse')
    indvar_subs(q) = mean(indvar);
    depvar_subs(q) = mean(depvar);
    else
        indvar_subs(:,q) = indvar;
        depvar_subs(q) = mean(depvar);
    end
    
    
    
    
    currentData = [];
end


if EasyParse(varargin,'RemoveOutliers','on')
    rmoutliers = [find(isoutlier(indvar_subs)) find(isoutlier(depvar_subs))];
    
    indvar_subs(rmoutliers) = [];
    depvar_subs(rmoutliers) = [];
end

if ~strcmpi(indvar_name,'itc_timecourse')
if ~strcmpi(indvar_name,'phase')
    [r,p] = corr(indvar_subs',depvar_subs','Type','Spearman');
    stats.rho = r;
    stats.p = p;
elseif strcmpi(depvar_name,'phase')
    [r,p] = circ_corrcl(depvar_subs,indvar_subs);
    stats.rho = r;
    stats.p = p;
else
    [r,p] = circ_corrcl(indvar_subs,depvar_subs);
    stats.rho = r;
    stats.p = p;
end

try
    if ~EasyParse(varargin,'Plot','off')
        if ~strcmpi(indvar_name,'phase') && ~strcmpi(depvar_name,'phase')
            nicecorrplot(indvar_subs,depvar_subs,{indvar_name depvar_name})
        else
            polarscatter(NormOntoRange(indvar_subs,[0 1]),depvar_subs,'filled')
        end
    end
catch
    warning('Plotting failed')
end
else
    stats = NaN;
end



