function [OsciFrac] = FracPower_EEG_wrapper(spec,frange,normalize)

OsciFrac = zeros(1,size(spec.osci,2));

if nargin < 3
   normalize = 0; 
end

if nargin < 2
   frange = [0.5 50]; 
end

disp(' ')
disp('Computing oscillatory/fractal ratio...')

freqindex = intersect(find(spec.freq > frange(1)),find(spec.freq < frange(2)));
freqs = spec.freq(freqindex);
    
for c = 1:size(spec.frac,2)
    fprintf([num2str(c) ' ']);
    %osci = sgolayfilt(spec.osci(:,c),5,1501);
    %osci = osci(freqindex);
    %frac = sgolayfilt(spec.frac(:,c),5,1501);
    frac = spec.mixd(freqindex,c);
    
    OsciFrac(c) = trapz(freqs,frac);
    if normalize
        OsciFrac = NormOntoRange(OsciFrac,[0 1]);
    end
end