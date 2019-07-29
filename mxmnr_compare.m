function [mdl1,mdl2,comparison] = mxmnr_compare(data,frmla1,frmla2)
% compares two logistic mixed-effects models using brms and loo in R
% data should be a table in this case with all relevant data included
% frmla1 and frmla2 are formulas for the two models


% if ~istable(X)
% xnames = cellstr(num2str([1:size(X,2)]'))';
% xnames = cellcat('X',xnames,'',0);
% X = array2table(X,'VariableNames',xnames);
% end
% 
% if ~istable(Y) && ~isempty(Y)
%     ynames = cellstr(num2str([1:size(Y,2)]'))';
%     ynames = cellcat('Y',xnames,'',0); 
%     Y = array2table(Y,'VariableNames',ynames);
% end
% 
% if ~isempty(Y)
%     datatbl = cat(2,X,Y);
% else
%     datatbl = X;
% end

currdir = pwd;

data{:,1:end-1} = zscore(data{:,1:end-1},[],1);

writetable(data,fullfile(currdir,'datatbl.csv'))

path = which('mxmnr_compare');

path = erase(path,'mxmnr_compare.m');

funcname = 'mxmnr_compare.R';
filein = 'datatbl.csv';
fileout = 'output';

[~,sysname] = system('hostname');

if contains(sysname,'Mac-the-Knife')
    setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
end

%system('ulimit -s 16384')

system(['R -e ''source("' fullfile(path,funcname) '"); mxmnr_compare("'...
    fullfile(currdir,filein) '",' frmla1 ',' frmla2 ',"' ... 
    fullfile(currdir,fileout) '")'''])

mdl1 = jsonread(fullfile(currdir,[fileout '_summary1.json']));
mdl2 = jsonread(fullfile(currdir,[fileout '_summary2.json']));
%comparison = jsonread(fullfile(currdir,[fileout '_cmpare.json']));

system(['rm ' filein]);
system(['rm ' fileout '_summary1.json']);
system(['rm ' fileout '_summary2.json']);
system(['rm ' fileout '_cmpare.json']);





