function [irasatc,time] = IRASA_power_timecourse(datain,srate,bandpass,oscifrac,savename)

rmpath(genpath('/home/sorenwt/projects/def-gnorthof/sorenwt/fieldtrip-master'))

if size(datain,1) > size(datain,2)
    datain = datain';
end

%if nargin < 6
    winsize = round(srate*3);
    overlap = 0.95;
%end

%if nargin < 5
%    specs = struct;
%end

if ~isfile(savename)
    calcSpecs = 1;
    specs = struct;
else
    calcSpecs = 0;
	load(savename)
end

if isempty(bandpass)
    bandpass = [1.3 85]; %hard-coded for HCP data
    %	oscifrac = 'frac';
end
disp(bandpass)


disp('Getting IRASA envelope...')
for c = 1:size(datain,1)
    fprintf([num2str(c) ' '])
if calcSpecs
        [datamat,time] = WindowToMatrix(datain(c,:),winsize,overlap);
        tmp = amri_sig_fractal(datamat,srate,'hset',[1.1:0.05:1.95 2.05:0.05:2.9]);
        if c == 1
            specs = tmp;
        else
            specs(c) = tmp;
        end
    end
    bpindex = intersect(find(specs(c).freq > bandpass(1)),find(specs(c).freq < bandpass(2)));
    irasatc(c,:) = trapz(specs(c).(oscifrac)(bpindex,:));
end

if calcSpecs
save(savename,'specs','time','-v7.3')
end

%all hardcoded for HCP data
dsampfact = 5;
irasatc = downsample(irasatc',dsampfact)';
time = downsample(time,dsampfact);

time = time/srate;

addpath('/home/sorenwt/projects/def-gnorthof/sorenwt/fieldtrip-master')
ft_defaults
%addpath(genpath(toolboxdir('signal')))
