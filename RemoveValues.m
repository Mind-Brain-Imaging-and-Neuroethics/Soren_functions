function [outthing] = RemoveValues(values,dim,inthing)
%values should be a vector of 0 and 1 containing the indices of the things
%you want to remove


indx = [];
for cc = 1:length(values)
    if values(cc)
        indx = [indx cc];
    end
end




outthing = inthing
if dim == 1
    outthing(indx,:) = [];
elseif dim == 2
    outthing(:,indx) = [];
end

