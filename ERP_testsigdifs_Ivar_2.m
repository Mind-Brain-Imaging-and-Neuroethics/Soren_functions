function [erp1,erp2,erpsub,times,sigdifmat] = ERP_testsigdifs_Ivar(conditions,varargin)
%currently you can only do pairwise comparisons

eeglab
directory_name =uigetdir;
cd(directory_name);

files = dir(fullfile(directory_name, '*.set'));

fileindex = find(~[files.isdir]);

filenames = extractfield(files,'name');

for c = 1:length(filenames)
   conditionnames{c} = extractBetween(filenames{c},'_','.set');
end

for c = 1:length(filenames)
    subnames{c} = extractBefore(filenames{c},'_');
end

[indices,allindices] = Ivar_findConds(conditions,conditionnames);

if EasyParse(varargin,'TrialBased','true')
    defEEG = pop_loadset('filename',filenames(1),'filepath',directory_name);
    for c = 1:2
        CurrEEG = eeg_emptyset;
        CurrEEG.data = [];
        CurrEEG.setname = ['Condition ' num2str(c)];
        CurrEEG.filepath = directory_name;
        CurrEEG.filename = ['Condition ' num2str(c) '.set'];
        CurrEEG.chanlocs = defEEG.chanlocs([9 10 52 24]);
        CurrEEG.srate = 500;
        CurrEEG.xmin = defEEG.xmin;
        CurrEEG.xmax = defEEG.xmax;
        for cc = 1:length(indices{c})
            EEG = pop_loadset('filename',filenames{indices{c}(cc)},'filepath',directory_name);
            if strcmpi(EEG.chanlocs(1).labels,'Fp1')
                CurrEEG.data = cat(3,CurrEEG.data,EEG.data([9 10 52 24],:,:)); %normal EEG cap
            else
                CurrEEG.data = cat(3,CurrEEG.data,EEG.data([19 20 40 51],:,:)); %fnirs EEG cap
            end
        end
        CurrEEG = eeg_checkset(CurrEEG);
        ALLEEG(c) = CurrEEG;
    end
    
    if ~CheckInput(varargin,'electrode')
        [erp1,erp2,erpsub,times,~] = pop_comperp(ALLEEG,1,1,2,'addavg','on','subavg','on');
    else
        [erp1,erp2,erpsub,times,~] = pop_comperp(ALLEEG,1,1,2,'chans',EasyParse(varargin,'electrode'),'addavg','on','subavg','on');
    end

    for c = 1:size(ALLEEG(1).data,2)
       for cc = 1:4
          sigdifmat(cc,c) = ranksum(squeeze(ALLEEG(1).data(cc,c,:)),squeeze(ALLEEG(1).data(cc,c,:)));
       end
    end

else 
    for c = 1:length(allindices)
        EEG = pop_loadset( 'filename', filenames{allindices(c)}, 'filepath', directory_name);
        ALLEEG(c) = EEG;
    end

    deletedsubs = {};
    removesubs = [];
    if iscell(indices{2})
        for c = 1:length(indices{1})
            sub1 = subnames{indices{2}{1}(c)};
            findsub = find(strcmpi(sub1,subnames(indices{2}{2})));
            set1 = find(allindices == indices{2}{1}(c));
            
            if ~isempty(findsub)
                set2 = find(allindices == indices{2}{2}(findsub));
                ALLEEG(set1).data = cat(3,ALLEEG(set1).data,ALLEEG(set2).data);
                ALLEEG(set1) = eeg_checkset(ALLEEG(set1));
                removesubs = [removesubs set2];
            else
                removesubs = [removesubs set1];
                deletedsubs = [deletedsubs {sub1}];
            end
        end
    end
    
    if ~isempty(deletedsubs)
        for c = 1:length(deletedsubs)
            for cc = 1:length(indices{1})
                if strcmpi(deletedsubs{c},subnames{indices{1}(cc)})
                    removesubs = [removesubs find(allindices == indices{1}(cc))]
                end
            end
        end
    end
    
    ALLEEG(removesubs) = [];
    
    if iscell(indices{2})
        conditions{2} = conditions{2}{1};
    end
    
    alleeg_conds = extractfield(ALLEEG,'filename');
    for c = 1:length(ALLEEG)
        alleeg_conds{c} = cellstr(extractBetween(alleeg_conds{c},'_','.set'));
    end
    alleeg_indices = Ivar_findConds(conditions,alleeg_conds)
    
    for c = 1:length(ALLEEG)
        if strcmpi(ALLEEG(c).chanlocs(1).labels,'Fp1')
            ALLEEG(c).data = ALLEEG(c).data([9 10 52 24],:,:); %normal EEG cap
            ALLEEG(c).chanlocs = ALLEEG(c).chanlocs([9 10 52 24]);
        else
            ALLEEG(c).data = ALLEEG(c).data([19 20 40 51],:,:); %fnirs EEG cap
            ALLEEG(c).chanlocs = ALLEEG(c).chanlocs([19 20 40 51]);
        end
        ALLEEG(c) = eeg_checkset(ALLEEG(c));
    end
    
    if ~CheckInput(varargin,'electrode')
        [erp1,erp2,erpsub,times,sigdifmat] = pop_comperp(ALLEEG,1,alleeg_indices{1},alleeg_indices{2},'alpha',0.05,'addavg','on','subavg','on');
    else
        [erp1,erp2,erpsub,times,sigdifmat] = pop_comperp(ALLEEG,1,alleeg_indices{1},alleeg_indices{2},'alpha',0.05,'chans',EasyParse(varargin,'electrode'),'addavg','on','subavg','on');
    end
end



if CheckInput(varargin,'Mcorrect')
    if EasyParse(varargin,'Mcorrect','fdr')
        if CheckInput(varargin,'electrode')
            sigdifmat = sigdifmat(EasyParse(varargin,'electrode'),1350:1650);
        else
            sigdifmat = sigdifmat(:,1350:1650);
        end
        
        sigdifmatdims = size(sigdifmat);
        sigdifmat = reshape(mafdr(reshape(sigdifmat,[],1),'BHFDR',true),sigdifmatdims(1),sigdifmatdims(2));
    end
end