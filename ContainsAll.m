function [torf] = ContainsAll(teststr,parts)

torf = 1;

for c = 1:length(parts)
   if ~contains(teststr,parts{c})
       torf = 0;
   end
end