function [dataout] = signrank_mat(datain1,datain2,dim)

if dim == 1
    datain1 = datain1';
    datain2 = datain2';
end

for c = 1:size(datain1,1)
if ~isempty(find(~isnan(datain1(c,:)))) && ~isempty(find(~isnan(datain2(c,:)))) 
  dataout(c) = signrank(datain1(c,:),datain2(c,:));
else
dataout(c) = NaN;
end
end
