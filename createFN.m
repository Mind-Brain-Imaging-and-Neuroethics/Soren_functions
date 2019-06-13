function signal = createFN(ple,siglength,method)

if nargin > 2 && (strcmpi(method,'r')) && siglength <= 60000
setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')

system(['R -e ' char(39) 'source("~/Documents/MATLAB/Functions/createFN.R"); createFN(' num2str(ple) ',' num2str(siglength) ')' char(39)])

tbl = readtable('/Users/Soren/output.csv');

signal = tbl.x;

!rm /Users/Soren/output.csv
else
    if ple < 1
        signal = MakeFGN(ple,siglength);
    else
        signal = wfbm((ple-1)/2,siglength);
    end
end