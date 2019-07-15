function [med,modtest] = ModeratedMediation(X,Y,M,W)

datatbl = array2table([vert(X),vert(Y),vert(M),vert(W)],'VariableNames',{'X','Y','M','W'});

writetable(datatbl,'/Users/Soren/datatbl.csv');

setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')

system('R -e ''source("~/Documents/MATLAB/Functions/ModMed.R"); ModMed("/Users/Soren/datatbl.csv","/Users/Soren/output")''')

med = jsonread('/Users/Soren/output_model.json');
modtest = jsonread('/Users/Soren/output_modtest.json');