function [output] = NormOntoRange(data,range,dimn)
%works along columns by default

if nargin < 3
    dimn = 1;
end

if dimn == 2
   data = data'; 
end

if size(data,1) == 1
   data = data'; 
end

output = zeros(size(data));

maxV = max(data,[],1);
minV = min(data,[],1);

for q = 1:size(data,2)
    for qq = 1:size(data,1)
        output(qq,q) = (range(2)-range(1))/(maxV(q)-minV(q))*(data(qq,q)-maxV(q))+range(2);
    end
end
