%function EpochGroups(EEG,subject,subID,boundaries,trialListenOrder)
%subject should be just one file in subjects, the one corresponding to the
%subject ID

%boundaries = GroupTrials(EEG)
%put this back in when you fix GroupTrials

groups = GetGroupBoundaries(subject(31:45,:));

for c = 1:80
    for cc = 1:14
        if isempty(find(subject(31:45,c) == 0)) %skip over if the person didn't group any tones
            if groups(cc,find(subject(50,:) == trialListenOrder(c)))
                if subject(49,find(subject(50,:) == trialListenOrder(c))) == 0
                    if subject(48,find(subject(50,:) == trialListenOrder(c))) == 0.7
                        EEG.event(boundaries(c)+cc*2-1).type = 'G1on_fast_05';
                    else
                        EEG.event(boundaries(c)+cc*2-1).type = 'G1on_fast_15';
                    end
                else
                    if subject(48,find(subject(50,:) == trialListenOrder(c))) == 0.7
                        EEG.event(boundaries(c)+cc*2-1).type = 'G1on_slow_05';
                    else
                        EEG.event(boundaries(c)+cc*2-1).type = 'G1on_slow_15';
                    end
                end
            elseif cc > 1 && groups(cc-1,trialListenOrder(c))
                if subject(49,find(subject(50,:) == trialListenOrder(c))) == 0
                    if subject(48,find(subject(50,:) == trialListenOrder(c))) == 0.7
                        EEG.event(boundaries(c)+cc*2-1).type = 'G2on_fast_05';
                    else
                        EEG.event(boundaries(c)+cc*2-1).type = 'G2on_fast_15';
                    end
                else
                    if subject(48,find(subject(50,:) == trialListenOrder(c))) == 0.7
                        EEG.event(boundaries(c)+cc*2-1).type = 'G2on_slow_05';
                    else
                        EEG.event(boundaries(c)+cc*2-1).type = 'G2on_slow_15';
                    end
                end
            else
                if subject(49,find(subject(50,:) == trialListenOrder(c))) == 0
                    if subject(48,find(subject(50,:) == trialListenOrder(c))) == 0.7
                        EEG.event(boundaries(c)+cc*2-1).type = 'NGon_fast_05';
                    else
                        EEG.event(boundaries(c)+cc*2-1).type = 'NGon_fast_15';
                    end
                else
                    if subject(48,find(subject(50,:) == trialListenOrder(c))) == 0.7
                        EEG.event(boundaries(c)+cc*2-1).type = 'NGon_slow_05';
                    else
                        EEG.event(boundaries(c)+cc*2-1).type = 'NGon_slow_15';
                    end
                end
            end
        end
    end
end

% rawconditions = {'G1on_fast_05','G1on_fast_15','G1on_slow_05','G1on_slow_15','G2on_fast_05','G2on_fast_15','G2on_slow_05','G2on_slow_15','NGon_fast_05','NGon_fast_15','NGon_slow_05','NGon_slow_15'};
% 
% conditions = [rawconditions {rawconditions(1:4)} {rawconditions(5:8)} {rawconditions(9:12)} {rawconditions([1 2 5 6 9 10])} {rawconditions([3 4 7 8 11 12])} {rawconditions([1 3 5 7 9 11])} {rawconditions([2 4 6 8 10 12])} ...
%     {rawconditions(1:2)} {rawconditions(3:4)} {rawconditions(5:6)} {rawconditions(7:8)} {rawconditions(9:10)} {rawconditions(11:12)} ...
%     {rawconditions([1 3])} {rawconditions([2 4])} {rawconditions([5 7])} {rawconditions([6 8])} {rawconditions([9 11])} {rawconditions([10 12])} ...
%     {rawconditions([1 5 9])} {rawconditions([2 6 10])} {rawconditions([3 7 11])} {rawconditions([4 8 12])}];
% 
% conditionNames = [rawconditions {'G1','G2','NG','Fast','Slow','cond_05','cond_15','G1on_fast','G1on_slow','G2on_fast','G2on_slow','NGon_fast','NGon_slow','G1on_05','G1on_15','G2on_05','G2on_15','NGon_05','NGon_15','Fast_05','Fast_15','Slow_05','Slow_15'}]; 
% 
% for c = 1:length(conditions)
%     if c <= 12
%         OUTEEG = pop_epoch(EEG,conditions(c),[-3 3]);
%     else
%         OUTEEG = pop_epoch(EEG,conditions{c},[-3 3]);
%     end
%     
%    OUTEEG = eeg_checkset(OUTEEG);
%    OUTEEG = pop_saveset(OUTEEG,'filename',[num2str(subID) '_' conditionNames{c}],'filepath','/Users/Soren/Desktop/SorenMusic/EEGData/EpochedReal/');
% end
% 


