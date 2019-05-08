function SelectTrials(EEG,subID,allself)
%EEG is the original data set for a given subject
%allself is the table of all the behavioural data


%currentSelf = allself(find(allself.id == subID),:);

%conditions are maincond(delay),matching,shapecat,and accuracy

%newsets = cell(4,2,3,3);

startpos = (subID-201)*400;

indx = [];
for c = 1:length(EEG.event)
   if strcmpi(EEG.event(c).type,'S  5')
      indx = [indx c];
   end
end

if length(indx) ~= 400
   disp('Buffer Overflow error: proceeding to next subject')
   return;
end

for c = 1:400
    EEG.event(indx(c)).type = [num2str(allself.maincond(startpos+c)) num2str(allself.matching(startpos+c)+1) num2str(allself.shape_cat(startpos+c)) num2str(allself.accuracy(startpos+c)+2)];
end

for c = 1:4
    for cc = 1:2
        for ccc = 1:3
            for cccc = 1:3
                try
                    newEEG = pop_selectevent(EEG,'type',[num2str(c) num2str(cc) num2str(ccc) num2str(cccc)],'deleteevents','off','deleteepochs','on','invertepochs','off');
                    newEEG.setname = [num2str(subID) '_' num2str(c) num2str(cc) num2str(ccc) num2str(cccc)];
                    newEEG = eeg_checkset(newEEG);
                    pop_saveset(newEEG,'filename',[newEEG.setname '.set'],'filepath','/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/raw_eeg/Self_conds');
                catch
                    disp(lasterror)
                end
            end
        end
    end
end


% for c = 1:size(EEG.data,3)
%    newsets{allself.maincond(startpos+c),allself.matching(startpos+c)+1,allself.shapecat(startpos+c),allself.accuracy(startpos+c)+2} = ...
%        cat(3,newsets{allself.maincond(startpos+c),allself.matching(startpos+c)+1,allself.shapecat(startpos+c),allself.accuracy(startpos+c)+2},EEG.data(:,:,c))
% end
