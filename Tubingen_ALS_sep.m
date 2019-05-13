function [als lis clis] = Tubingen_ALS_sep(outputs)

lissub = 'SELFA3';

clissub = {'SELFA2' 'SELFA8'};

allsubs = outputs.sub;
allsubs = extractBefore(allsubs,'S001');

lisindx = find(strcmpi(lissub,allsubs));

clisindx = Subset_index(allsubs,clissub);

lis = outputs;
lis.data = lis.data(lisindx,:,:);
lis.sub = lissub;

clis = outputs;
clis.data = clis.data(clisindx,:,:);
clis.sub = clissub;

rmindx = [lisindx clisindx];
alsindx = 1:length(outputs.sub);
alsindx = alsindx(find(~ismember(alsindx,rmindx)));
als = outputs;
als.data = als.data(alsindx,:,:);
als.sub = als.sub(alsindx);