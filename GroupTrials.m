function [boundaries,listenbounds] = GroupTrials(EEG)

boundaries = [];
ind = 1;

for c = 2:length(EEG.event)
    %if strcmpi(EEG.event(c).type,'S  1')
    %if EEG.event(c).latency - EEG.event(c-1).latency > EEG.srate && strcmpi(EEG.event(c).type,'S  1')
    %if ((EEG.event(c).latency - EEG.event(c-1).latency > EEG.srate) | (~strcmpi(EEG.event(c-1).type,'S  3'))) && strcmpi(EEG.event(c).type,'S  4')
    if strcmpi(EEG.event(c-1).type,'S  2') && strcmpi(EEG.event(c).type,'S  4')
        boundaries(ind) = c;
        ind = ind+1;
    end
end


%boundaries = boundaries - 1;

listenbounds = 0;

for c = 1:length(boundaries)-1
   if boundaries(c+1) - boundaries(c) < 50
       listenbounds = boundaries(c)
       break;
   end
end



