function windows = getWindows(time,winsize,toi,data)
% outputs a sliding window decomposition of optional 2D matrix data with
% the size specified by winsize and the times specified by toi. If data is
% not entered, outputs indices of time. Winsize should be specified in
% terms of units of time

%tindex = FindClosest(time,toi);

for c = 1:length(toi)
%try
    windows{c} = FindClosest(time,toi(c)-winsize*0.5):FindClosest(time,toi(c)+winsize*0.5);
% catch
%         test = length(FindClosest(time,toi(c)-winsize*0.5):FindClosest(time,toi(c)+winsize*0.5));
%         
%         if length(windows(c-1,:)) > test
%            windows(c,:) = FindClosest(time,toi(c)-winsize*0.5):FindClosest(time,toi(c)+winsize*0.5,2);
%         end
% end
end

if nargin > 3
   for c = 1:length(toi)
      windows{c} = data(:,windows{c}); 
   end
end