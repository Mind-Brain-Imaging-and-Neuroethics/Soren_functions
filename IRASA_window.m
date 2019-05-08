function [specout] = IRASA_window(datain,srate,varargin)
%datain should be in the form channels x time points


if ~CheckInput(varargin,'Winsize')
winsize = round(10*srate);
else
    winsize = EasyParse(varargin,'Winsize');
end

if ~CheckInput(varargin,'Overlap')
    overlap = 0;
else
    overlap = EasyParse(varargin,'Overlap');
end

if ~CheckInput(varargin,'hset')
    hset = [1.1:0.05:1.95 2.05:0.05:2.9];
else
    hset = EasyParse(varargin,'hset');
end

for c = 1:size(datain,1)
   specs(c) = amri_sig_fractal(WindowToMatrix(datain(c,:),winsize,overlap),srate,'hset',hset);
end

for c = 1:length(specs)
   specs(c).mixd = nan_mean(specs(c).mixd,2);
   specs(c).frac = nan_mean(specs(c).frac,2);
   specs(c).osci = nan_mean(specs(c).osci,2);
end

specout = specs(1);
specout.mixd = [];
specout.frac = [];
specout.osci = [];
for c = 1:length(specs)
   specout.mixd = cat(2,specout.mixd,specs(c).mixd);
   specout.frac = cat(2,specout.frac,specs(c).frac);
   specout.osci = cat(2,specout.osci,specs(c).osci);
end
end

