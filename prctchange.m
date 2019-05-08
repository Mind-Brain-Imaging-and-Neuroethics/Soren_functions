function [output] = prctchange(inmat)
%works along rows

output = inmat;

for c  = 1:size(inmat,1)
   output(c,:) = ((inmat(c,:)-mean(inmat(c,:)))*100)./mean(inmat(c,:));
end