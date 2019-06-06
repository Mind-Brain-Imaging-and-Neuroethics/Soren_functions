function [realsig] = Make_signal_to_ampenv(ampenv,period,offset)

if nargin < 3
   offset = 0; 
end

signal = linspace(offset,2*pi+offset,period+1);
signal(end) = [];

signal = repmat(signal,1,floor(length(ampenv)/length(signal)));
signal = [signal signal(1:(length(ampenv)-length(signal)))];

for c = 1:length(signal)
   [x,y] = pol2cart(signal(c),ampenv(c));
   realsig(c) = x;
end
