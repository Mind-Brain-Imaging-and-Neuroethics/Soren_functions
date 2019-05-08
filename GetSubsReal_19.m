function [outmat] = GetSubsReal_19(inmat,meandata)

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/trialFlags_ivar_2.mat')
trialFlags = rmfield(trialFlags,'sub214');
trialFlags = rmfield(trialFlags,'sub224');
fields = fieldnames(trialFlags);

for c = 1:length(fields)
    outmat(c,:) = inmat(find(meandata.subid == str2num(erase(fields{c},'sub'))),:);
end