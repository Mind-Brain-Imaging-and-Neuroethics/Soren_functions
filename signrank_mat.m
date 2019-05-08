function [dataout] = signrank_mat(datain1,datain2,dim)

if dim == 1
    datain1 = datain1';
    datain2 = datain2';
end

for c = 1:size(datain1,1)
   dataout(c) = signrank(datain1(c,:),datain2(c,:));
end