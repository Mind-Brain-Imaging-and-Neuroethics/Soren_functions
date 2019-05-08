function [pdifs,rawdifs] = LZC_comparecond_Ivar(LZCbefore,LZCafter,conditions,electrode,varargin)

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/trialFlags_ivar.mat')

fields = fieldnames(trialFlags);

data = cell(length(conditions),2);


for i = 1:length(fields)
    [~,~,flatindices] = Ivar_findConds_2(conditions,trialFlags.(fields{i}));
    for c = 1:length(conditions)
        data{c,1} = [data{c,1} LZCbefore.(fields{i})(electrode,flatindices{c})];
        data{c,2} = [data{c,2} LZCafter.(fields{i})(electrode,flatindices{c})];
    end
end

if ~EasyParse(varargin,'Mean','on')
    for c = 1:length(electrode)
        pdifs(c) = ranksum(data{1,1}(electrode(c),:),data{1,2}(electrode(c),:));
        rawdifs(c) = mean(data{1,1}(electrode(c),:)) - mean(data{1,2}(electrode(c),:));
    end
else
    pdifs = ranksum(mean(data{1,1},1),mean(data{1,2},1));
    rawdifs = mean(mean(data{1,1})) - mean(mean(data{1,2}));
end

if EasyParse(varargin,'Plot','on')
   boxplot(horzcat(mean(data{1,1},1)',mean(data{1,2},1)'),'Labels',{'Prestimulus LZC','Poststimulus LZC'});
   title(['LZC difference in condition ' conditions{1}])
end
