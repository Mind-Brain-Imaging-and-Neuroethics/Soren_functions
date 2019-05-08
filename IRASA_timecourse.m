function [irasatc,time,specs] = IRASA_power_timecourse(datain,srate,bandpass,oscifrac,winsize,overlap)

if size(datain,1) > size(datain,2)
   datain = datain'; 
end

if nargin < 3
   winsize = srate*3;
   overlap = srate*0.125;
end

for c = 1:length(datain)
    [datamat,time] = WindowToMatrix(datain(c,:),winsize,overlap)
    specs(c) = amri_sig_fractal(datamat,srate,'hset',1.1:0.05:2.9);
    bpindex = intersect(find(specs(c).freq > bandpass(1)),find(specs(c).freq < bandpass(2)));
    irasatc(c,:) = trapz(specs(c).(oscifrac)(bpindex,:));
end

time = time/srate;