function Topoplot_Ivar(datamat)
%datamat should be in format rows = subjects, columns = channels

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/chanlocs.mat')
load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/trialFlags_ivar_2.mat')
load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/meandata_2.mat')

trialFlags = rmfield(trialFlags,'sub214');
trialFlags = rmfield(trialFlags,'sub224');
fields = fieldnames(trialFlags);

fnirsmat = [];
normalmat = [];
for i = 1:length(fields)
    if meandata.fnirscap(find(meandata.subid == str2num(erase(fields{i},'sub'))))
       fnirsmat(i,:) = datamat(i,:);
    else
        normalmat(i,:) = datamat(i,:);
    end
end

subplot(1,2,1)
topoplot(mean(normalmat,1),normal_chanlocs)
subplot(1,2,2)
topoplot(mean(fnirsmat,1),fnirs_chanlocs)

