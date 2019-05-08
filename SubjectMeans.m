function [outvect] = SubjectMeans(inputmatrix)

numsubs = size(inputmatrix,1);
numelectrodes = size(inputmatrix,2);

tmp = reshape(inputmatrix',[],1);

for cc = 1:31
    outvect(cc) = mean(inputmatrix(((cc-1)*numelectrodes+1):(cc*numelectrodes)));
end

