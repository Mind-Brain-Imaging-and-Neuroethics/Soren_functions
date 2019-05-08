function [pdifs,datattv,tpoints,ttv1,ttv2,pseudoTTV] = TTV_compareconds(ttvvals,electrodes,conditions,varargin)

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/trialFlags_ivar_2.mat')

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/meandata_2.mat');

trialFlags = rmfield(trialFlags,'sub214');
trialFlags = rmfield(trialFlags,'sub224');

fields = fieldnames(trialFlags);

data = cell(1,2) ;

if CheckInput(varargin,'SlidingWindow')
    winparams = EasyParse(varargin,'SlidingWindow');
end

%tpoints = NaN;

if CheckInput(varargin,'Limits')
    tl = EasyParse(varargin,'Limits');
    for i = 1:length(fields)
        ttvvals.(fields{i}) = ttvvals.(fields{i})(:,tl(1):tl(2),:);
    end
end

if CheckInput(varargin,'PrestimSplit')
    splitparams = EasyParse(varargin,'PrestimSplit');
    %     prestim = EasyParse(splitparams,'Prestim')
    %
    %     if CheckInput(splitparams,'SlidingWindow')
    %         [~,datattv,tpointsTTV] = TTV_compareconds(prestim,1:63,{'xxxxx','43211'},'Plot','off','MCorrect','off','SigDifs','off','SlidingWindow',EasyParse(splitparams,'SlidingWindow'));
    %     else
    %         [~,datattv,tpointsTTV] = TTV_compareconds(prestim,1:63,{'xxxxx','43211'},'Plot','off','MCorrect','off','SigDifs','off');
    %     end
    %
    
    trange = EasyParse(splitparams,'Trange');
    frange = EasyParse(splitparams,'Frange');
    indvar_name = EasyParse(splitparams,'Criterion');

    [~,indvar_record] = TrialCorr_Ivar(trange,frange,electrodes,indvar_name,'RT');
    
    for i = 1:length(fields)
        if strcmpi(indvar_name,'phase')
            if EasyParse(splitparams,'Phasemethod','cos')
                indvar_record.(fields{i}) = cos(indvar_record.(fields{i}));
                tmp = (indvar_record.(fields{i}) > 0)+1;
                [~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(fields{i}));
                mask = zeros(1,400);
                mask(flatindices{1}) = 1;
                tmp = tmp.*mask;
                tmp = cellstr(num2str(tmp));
                tmp = tmp';
            elseif EasyParse(splitparams,'Phasemethod','concentrated')
                warning('Using concentrated angle')
                indvar_record.(fields{i}) = wrapTo2Pi(indvar_record.(fields{i}));
                [~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(fields{i}));
                
                cmean = circ_mean(indvar_record.(fields{i})(flatindices{1}));
                cmean = wrapTo2Pi(cmean);
                
                concrange = [cmean-pi/2 cmean+pi/2];
                concrange = wrapTo2Pi(concrange);
                
                antirange = [cmean+pi/2 cmean+3*pi/2];
                antirange = wrapTo2Pi(antirange);
                
                tmp = zeros(1,400);
                if concrange(2) > concrange(1)
                    tmp(intersect(find(indvar_record.(fields{i}) > concrange(1)),find(indvar_record.(fields{i}) < concrange(2)))) = 1;
                else
                    tmp(intersect(find(indvar_record.(fields{i}) > concrange(1)),find(indvar_record.(fields{i}) < (2*pi)))) = 1;
                    tmp(intersect(find(indvar_record.(fields{i}) < concrange(2)),find(indvar_record.(fields{i}) > 0))) = 1;
                end
                
                if antirange(2) > antirange(1)
                    tmp(intersect(find(indvar_record.(fields{i}) > antirange(1)),find(indvar_record.(fields{i}) < antirange(2)))) = 2;
                else
                    tmp(intersect(find(indvar_record.(fields{i}) > antirange(1)),find(indvar_record.(fields{i}) < (2*pi)))) = 2;
                    tmp(intersect(find(indvar_record.(fields{i}) < antirange(2)),find(indvar_record.(fields{i}) > 0))) = 2;
                end
                
                mask = zeros(1,400);
                mask(flatindices{1}) = 1;
                tmp = tmp.*mask;
                
                tmp = tmp';
                tmp = cellstr(num2str(tmp));
                tmp = tmp';
            end
        else
            [~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(fields{i}));
            
            tmp = (indvar_record.(fields{i}) > median(indvar_record.(fields{i})(flatindices{1})))+1;
            tmp = tmp';
            mask = zeros(1,400);
            mask(flatindices{1}) = 1;
            tmp = tmp.*mask;
            tmp = tmp';
            tmp = cellstr(num2str(tmp));
            tmp = tmp';
            %                  tmp = cellstr(num2str((indvar_record.(fields{i}) > median(indvar_record.(fields{i})))+1));
            %                  tmp = tmp';
        end
        
        trialFlags.(fields{i}) = strcat(trialFlags.(fields{i}),tmp);
        
        
    end
    
    orig_conditions = conditions;
    conditions{1} = strcat(conditions{1},'1');
    conditions{2} = strcat(conditions{2},'2');
    %     datattv = datattv{1};
    %     datattv = datattv(:,trange,:);
    %     datattv = mean(mean(datattv,2),1);
    %     datattv = squeeze(datattv);
    
    
    
    %clear prestim;
else
    orig_conditions = conditions;
end


if CheckInput(varargin,'Pseudotrial')
    if ~strcmpi(EasyParse(varargin,'Pseudotrial'),'broadband')
        load(['/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTV' EasyParse(varargin,'Pseudotrial') '.mat'],'prestim')
%         if EasyParse(varargin,'Pseudotrial','alpha')
%             prestim = prestim_alpha;
%             clear prestim_alpha;
%         elseif EasyParse(varargin,'Pseudotrial','beta')
%             prestim = prestim_beta;
%             clear prestim_beta;
%         end
    else
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','prestim');
    end
    for i = 1:length(fields)
        prestim.(fields{i}) = prestim.(fields{i})(:,501:1000,:);
        if tl(2) - tl(1) < 499
            r = tl(2) - tl(1);
            prestim.(fields{i}) = prestim.(fields{i})(:,(500-r):500,:);
        end
    end
    
           %warning('Using raw conditions - if these have been altered due to a prestim split, errors may occur or results may be inaccurate')

    for i = 1:length(fields)
       [~,~,flatindices] = Ivar_findConds_2(orig_conditions,trialFlags.(fields{i}));
       pseudoTTV{1}(:,:,i)  = std(prestim.(fields{i})(:,:,flatindices{1}),0,3); 
       pseudoTTV{2}(:,:,i) = std(prestim.(fields{i})(:,:,flatindices{2}),0,3); 
    end
    
elseif CheckInput(varargin,'PseudotrialSplit')
    if strcmpi(EasyParse(varargin,'PseudotrialSplit'),'broadband')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVcalc.mat','prestim');
    elseif strcmpi(EasyParse(varargin,'PseudotrialSplit'),'prestim2')
        load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat','prestim2');
        prestim = prestim2;
        clear prestim2
    else
        load(['/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTV' EasyParse(varargin,'PseudotrialSplit') '.mat'],'prestim')
%         if EasyParse(varargin,'PseudotrialSplit','alpha')
%             prestim = prestim_alpha;
%             clear prestim_alpha;
%         elseif EasyParse(varargin,'PseudotrialSplit','beta')
%             prestim = prestim_beta;
%             clear prestim_beta;
%         end
    end
    
    
    for i = 1:length(fields)
        [~,~,flatindices] = Ivar_findConds_2(orig_conditions,trialFlags.(fields{i}));
        prestim.(fields{i}) = prestim.(fields{i})(:,501:1000,flatindices{1});
    end
    %     if tl(2) - tl(1) < 499
    %         r = tl(2) - tl(1);
    %         prestim.(fields{i}) = prestim.(fields{i})(:,(500-r):500,:);
    %     end
    if ismatrix(EasyParse(varargin,'Pseudocrit'))
        [~,ersprecord] = TrialCorr_Ivar([-1100 1000],EasyParse(varargin,'Pseudocrit'),electrodes,'raw_ersp','RT','Conditions',orig_conditions,'LongEpochs','true');
        disp('Fix conditions stuff later')
    end
    for i = 1:length(fields)
        if EasyParse(varargin,'Pseudocrit','amp')
            pseudocrit = hilbert(t3d(prestim.(fields{i})));
            pseudocrit = squeeze(mean(abs(pseudocrit(1:50,:,:)),1));
        elseif ismatrix(EasyParse(varargin,'Pseudocrit'))
            pseudocrit = ersprecord.(fields{i});
            for c = 1:63
               pc(c,:) = pseudocrit;
            end
            pseudocrit = pc;
            clear pc
            warning('Method will fail if using multiple electrodes')
        end
%         if ~CheckInput(varargin,'Pseudocrit')
%             for c = 1:63
%                 pseudoTTV{1}(:,:,i) = std(prestim.(fields{i})(c,51:end,find(pseudocrit < median(pseudocrit))),0,3);
%                 pseudoTTV{2}(:,:,i) = std(prestim.(fields{i})(c,51:end,find(pseudocrit > median(pseudocrit))),0,3);
%             end
%         else
            for c = 1:63
                pseudoTTV{1}(c,:,i) = std(prestim.(fields{i})(c,(501-tl(2)):end,find(pseudocrit(c,:) < median(pseudocrit(c,:)))),0,3);
                pseudoTTV{2}(c,:,i) = std(prestim.(fields{i})(c,(501-tl(2)):end,find(pseudocrit(c,:) > median(pseudocrit(c,:)))),0,3);
            end
        %end
        
    end
end

if CheckInput(varargin,'TrimSplit') && EasyParse(EasyParse(varargin,'TrimSplit'),'Criterion','amplitude')
    load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/TTVprestim2.mat')
end

for i = 1:length(fields)
    disp(i)
    %disp(['Processing subject ' extractAfter(fields{i},'sub')])
    [~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(fields{i}));
    
    if CheckInput(varargin,'TrimSplit')
        splitparams = EasyParse(varargin,'TrimSplit');
        
        if EasyParse(splitparams,'Criterion','amplitude')
            tprange = EasyParse(splitparams,'Tprange');
            %prestim2.(fields{i}) = prestim2.(fields{i})(:,:,flatindices{1});
            for c = 1:length(flatindices{1}) %requires that the conditions be the same - just a pruning thing
                tmp = abs(hilbert(squeeze(mean(prestim2.(fields{i})(:,:,flatindices{1}(c)),1))));
                tmp = tmp(tprange(1):tprange(2));
                tmp = mean(tmp);
                prestimvals(c) = tmp;
            end
            
            flatindices{2} = flatindices{1}(find(prestimvals > median(prestimvals)));
            flatindices{1} = flatindices{1}(find(prestimvals < median(prestimvals)));
            clear prestimvals
        end
    end
    
    minTrials = min(cellfun('length',flatindices)); %this is the number of trials for the condition with the fewest trials
    numTrials = cellfun('length',flatindices);
    minIndex = find(numTrials == minTrials);
    
    if minTrials > -1
        for c = 1:length(conditions)
            %data{c} = cat(3,data{c},ttvvals.(fields{i})(:,:,flatindices{c}));
            
            if CheckInput(varargin,'NormTrials') && length(flatindices{c}) > minTrials
                for q = 1:EasyParse(varargin,'NormTrials')
                    bootIndices = randsample(flatindices{c},minTrials);
                    %bootstrapTTV(:,:,q) = Real_nanstd(ttvvals.(fields{i})(:,:,bootIndices),3);
                    bootstrapTTV(:,:,q,i) = std(ttvvals.(fields{i})(:,:,bootIndices),0,3);
                end
                datattv{c}(:,:,i) = mean(squeeze(bootstrapTTV(:,:,:,i)),3);
            else
                %datattv{c}(:,:,i) = Real_nanstd(ttvvals.(fields{i})(:,:,flatindices{c}),3);
                datattv{c}(:,:,i) = std(ttvvals.(fields{i})(:,:,flatindices{c}),0,3);
            end
            if CheckInput(varargin,'NormRange')
                datattv{c}(:,:,i) = ((datattv{c}(:,:,i)-mean(datattv{c}(:,EasyParse(varargin,'NormRange'),i),2))*100)./mean(datattv{c}(:,EasyParse(varargin,'NormRange'),i),2);
            elseif CheckInput(varargin,'Pseudotrial')
                warning('Dividing by pseudotrial mean')
                datattv{c}(:,:,i) = ((datattv{c}(:,:,i)-pseudoTTV{c}(:,:,i))*100)./(pseudoTTV{c}(:,:,i));
            elseif CheckInput(varargin,'PseudotrialSplit')
                warning('Dividing by pseudotrial mean')
                datattv{c}(:,:,i) = ((datattv{c}(:,:,i)-mean(pseudoTTV{c}(:,:,i),2))*100)./mean(pseudoTTV{c}(:,:,i),2);
            end
        end
    else
        datattv{1}(:,:,i) = NaN(size(ttvvals.(fields{i}),1),size(ttvvals.(fields{i}),2));
        datattv{2}(:,:,i) = NaN(size(ttvvals.(fields{i}),1),size(ttvvals.(fields{i}),2));
        warning('Too few trials - skipping subject')
    end
    
    try
        fnirscap(i) = meandata.fnirscap(find(meandata.subid == str2num(extractAfter(fields{i},'sub'))));
    catch
        disp(fields{i})
    end
end


%     for i = 1:length(fields)
%        for c = 1:length(conditions)
%           for cc = 1:size(datattv{c},1)
%               newdatattv{c}(cc,:,i) = datattv{c}(cc,:,i) - Real_nanstd(prestim2.(fields{i})(cc,:,:),3)
%           end
%        end
%     end
% end

if CheckInput(varargin,'SlidingWindow')
    for i = 1:length(fields)
        for c = 1:length(conditions)
            for cc = 1:size(datattv{c},1)
                [newdatattv{c}(cc,:,i) tpoints] = SlidingWindow(@nan_std,datattv{c}(cc,:,i),winparams(1),winparams(2));
            end
        end
    end
    
    datattv = newdatattv;
end

if CheckInput(varargin,'SlidingWindow')
    varargin = [varargin {'Xaxis',tpoints*2}];
end

pdifs = 0;

elecindex.fnirs = [13 14 19 20 40 51 45 46];

elecindex.normal = [37 38 9 10 52 24 56 57];

elecnames.fnirs = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};

elecnames.normal = {'F1','F2','Fc1','Fc2','Cpz','Pz','P1','P2'};

if iscell(electrodes)
    newdatattv = [];
    for i = 1:length(fields)
        if fnirscap(i)
            newdatattv{1}(:,:,i) = datattv{1}(elecindex.fnirs(find(contains(elecnames.fnirs,electrodes))),:,i);
            newdatattv{2}(:,:,i) = datattv{2}(elecindex.fnirs(find(contains(elecnames.fnirs,electrodes))),:,i);
            if CheckInput(varargin,'Pseudotrial') || CheckInput(varargin,'PseudotrialSplit')
                newpseudottv{1}(:,:,i) = pseudoTTV{1}(elecindex.fnirs(find(contains(elecnames.fnirs,electrodes))),:,i);
                newpseudottv{2}(:,:,i) = pseudoTTV{2}(elecindex.fnirs(find(contains(elecnames.fnirs,electrodes))),:,i);
            end
        else
            newdatattv{1}(:,:,i) = datattv{1}(elecindex.normal(find(contains(elecnames.normal,electrodes))),:,i);
            newdatattv{2}(:,:,i) = datattv{2}(elecindex.normal(find(contains(elecnames.normal,electrodes))),:,i);
            if CheckInput(varargin,'Pseudotrial') || CheckInput(varargin,'PseudotrialSplit')
                newpseudottv{1}(:,:,i) = pseudoTTV{1}(elecindex.normal(find(contains(elecnames.normal,electrodes))),:,i);
                newpseudottv{2}(:,:,i) = pseudoTTV{2}(elecindex.normal(find(contains(elecnames.normal,electrodes))),:,i);
            end
        end
    end
    
    datattv = newdatattv;
    if CheckInput(varargin,'Pseudotrial') || CheckInput(varargin,'PseudotrialSplit')
       pseudoTTV = newpseudottv; 
    end
    electrodes = 1:length(electrodes);
end


if ~EasyParse(varargin,'SigDifs','off')
    if ~CheckInput(varargin,'NormTrials')
    for c = 1:size(datattv{1},2)
        pdifs(c) = signrank(squeeze(nan_mean(datattv{1}(electrodes,c,:),1)),squeeze(nan_mean(datattv{2}(electrodes,c,:),1)));
    end
    else
        for c = 1:size(datattv{minIndex},2)
           for cc = 1:EasyParse(varargin,'NormTrials')
              tmpDifs(c,cc) = signrank(squeeze(nan_mean(datattv{minIndex}(electrodes,c,:),1)),squeeze(nan_mean(bootstrapTTV(electrodes,c,cc,:),1)));
           end
           pdifs(c) = Fisher_combine(tmpDifs(c,:));
        end
    end
end

if ~EasyParse(varargin,'MCorrect','off')
    if sum(pdifs > 0.05) > 0 || sum(pdifs > 0.05) == length(pdifs)
        pdifs = mafdr(pdifs,'BHFDR',true);
    else
        pdifs = mafdr(pdifs,'BHFDR',true)
    end
end

ttv1 = nan_mean(nan_mean(datattv{1}(electrodes,:,:),1),3);
ttv2 = nan_mean(nan_mean(datattv{2}(electrodes,:,:),1),3);

xax = linspace(1,2*length(ttv1),length(ttv1));

if CheckInput(varargin,'Xaxis')
    xax = EasyParse(varargin,'Xaxis');
end

tpoints = xax;

if ~EasyParse(varargin,'Plot','off')
    if ~EasyParse(varargin,'Plot','all')
        figure
        if ~EasyParse(varargin,'PlotPseudo','true')
            
            plot(xax,ttv1,'b');
            
            hold on;
            plot(xax,ttv2,'r');
        else
            plot(xax,nan_mean(nan_mean(((datattv{1}(electrodes,:,:).*pseudoTTV{1}(electrodes,:,:)/100)+pseudoTTV{1}(electrodes,:,:)),1),3),'b');
            hold on
            plot(xax,nan_mean(nan_mean(((datattv{2}(electrodes,:,:).*pseudoTTV{2}(electrodes,:,:)/100)+pseudoTTV{2}(electrodes,:,:)),1),3),'r');            
            hold on
            plot(xax,nan_mean(nan_mean(pseudoTTV{1}(electrodes,:,:),1),3),'b--');
            hold on
            plot(xax,nan_mean(nan_mean(pseudoTTV{2}(electrodes,:,:),1),3),'r--');
        end
        
        patchindex = pdifs < 0.05;
        yl = ylim;
        if ~EasyParse(varargin,'SigDifs','off')
            patchstep = patchindex*(yl(2)-yl(1));
            patchstep = patchstep + yl(1);
            area(xax,patchstep,yl(1),'FaceAlpha',0.3,'LineStyle','none','FaceColor',[0.5 0.5 0.5])
            if CheckInput(varargin,'Legend')
                legend([EasyParse(varargin,'Legend') {'Significant Differences'}])
            else
                legend({'Condition 1','Condition 2','Significant Differences'})
            end
            ylim(yl)
        else
            if CheckInput(varargin,'Legend')
                legend([EasyParse(varargin,'Legend') {'Significant Differences'}])
            else
                legend({'Condition 1','Condition 2','Significant Differences'})
            end
        end
        xlabel('Time (ms)')
        if CheckInput(varargin,'NormRange') || CheckInput(varargin,'Pseudotrial') || CheckInput(varargin,'PseudotrialSplit')
            if CheckInput(varargin,'SlidingWindow')
                ylabel('SD of % change in TTV')
            else
                ylabel('% change in trial-to-trial SD')
            end
        else
            if CheckInput(varargin,'SlidingWindow')
                ylabel('SD of TTV')
            else
                ylabel('Trial-to-trial SD')
            end
        end
    else
        for q = 1:length(fields)
            figure
            
            plot(xax,nan_mean(datattv{1}(electrodes,:,q),1));
            
            hold on;
            plot(xax,nan_mean(datattv{2}(electrodes,:,q),1));
            %         patchindex = pdifs < 0.05;
            %         yl = ylim;
            %         patchstep = patchindex*(yl(2)-yl(1));
            %         patchstep = patchstep + yl(1);
            %         %    for c = 1:length(patchstep)
            %         %        patchstep2(((c-1)*2)+1) = patchstep(c);
            %         %        patchstep2(c*2) = patchstep(c);
            %         %    end
            %         area(xax,patchstep,yl(1),'FaceAlpha',0.3,'LineStyle','none','FaceColor',[0.5 0.5 0.5])
            if CheckInput(varargin,'Legend')
            legend([EasyParse(varargin,'Legend') {'Significant Differences'}])
            else
                legend({'Condition 1','Condition 2','Significant Differences'})
            end
            %        ylim(yl)
            xlabel('Time (ms)')
            ylabel('Trial-to-trial SD')
        end
    end
    
    %    for c = 1:length(ttv1)
    %        patch([t(patchindex) fliplr(t(patchindex))], [ones(size(ttv1(patchindex)))*yl(1) ones(size(ttv1(patchindex)))*yl(2)], [0.6 0.4 0.9], 'FaceAlpha',0.3, 'EdgeColor','none')
    %    end
end