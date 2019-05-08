function [output] = GetSubsReal(input,EEGM,meandata)

thingy = [];
for c = 1:height(EEGM)
   if mod(c,5) == 1
        thingy = [thingy; EEGM(c,:)];
   end
end

output = input;

indx = [];
for c = 1:44
   if thingy.restses(c) == 2
       indx = [indx c];
   end
end
% 
%output(indx,:) = [];
thingy(indx,:) = [];
% 
output = [thingy.subject output];

output = sortrows(output,1);


indx = [];
for c = 1:size(output,1)
   if isempty(find(meandata.subid == output(c,1)))
       indx = [indx c];
   end
end

output(indx,:) = [];

output = output(:,2:end);
