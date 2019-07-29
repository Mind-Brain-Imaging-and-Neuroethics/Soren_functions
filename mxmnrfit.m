function [mdl,summary] = CVtest(X,Y,frmla)
% performs a multinomial mixed-effects logistic regression on the data in
% Y, given the predictors in X
%
% X and Y can be tables with predictor names, or just matrices (in this case the
% names are input as X1,X2,etc)
% if Y is entered as empty, it is assumed that the response values are
% present in X


if ~istable(X)
xnames = cellstr(num2str([1:size(X,2)]'))';
xnames = cellcat('X',xnames,'',0);
X = array2table(X,'VariableNames',xnames);
end

if ~istable(Y) && ~isempty(Y)
    ynames = cellstr(num2str([1:size(Y,2)]'))';
    ynames = cellcat('Y',xnames,'',0); 
    Y = array2table(Y,'VariableNames',ynames);
end

if ~isempty(Y)
    datatbl = cat(2,X,Y);
else
    datatbl = X;
end


currdir = pwd;

writetable(datatbl,fullfile(currdir,'datatbl.csv'))

path = which('mxmnrfit');

path = erase(path,'mxmnrfit.m');

funcname = 'mxmnrfit.R';
filein = 'datatbl.csv';
fileout = 'output';

[~,sysname] = system('hostname');

if contains(sysname,'Mac-the-Knife')
    setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
end

system('ulimit -s 16384')

system(['R -e ''source("' fullfile(path,funcname) '"); mxmnrfit("'...
    fullfile(currdir,filein) '",' frmla ',"' ... 
    fullfile(currdir,fileout) '")'''])

mdl = jsonread(fullfile(currdir,[fileout '_model.json']));
summary = jsonread(fullfile(currdir,[fileout '_summary.json']));

system(['rm ' filein]);
system(['rm ' fileout '_model.json']);
system(['rm ' fileout '_summary.json']);





