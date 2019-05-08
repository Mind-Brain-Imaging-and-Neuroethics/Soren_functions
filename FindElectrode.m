function [elecindex] = FindElectrode(subNum,electrode,chanlocs,meandata,trialFlags)

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/chanlocs.mat')

load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/meandata_2.mat')

if subNum < 100 && nargin < 5
   load('/Volumes/SOREN_SSD/Stuff_Ivar_SSD/eeg/trialFlags_ivar_2.mat')
   fields = fieldnames(trialFlags);
   subNum = fields{subNum};
   subNum = erase(subNum,'sub');
   subNum = str2num(subNum);
end

fnirs = meandata.fnirscap(find(meandata.subid == subNum));

if fnirs
   elecs = extractfield(fnirs_chanlocs,'labels');
   for c = 1:63
      if strcmpi(electrode,elecs{c})
          elecindex = c;
          break
      end
   end
else
       elecs = extractfield(normal_chanlocs,'labels');
   for c = 1:63
      if strcmpi(electrode,elecs{c})
          elecindex = c;
          break
      end
   end
end
