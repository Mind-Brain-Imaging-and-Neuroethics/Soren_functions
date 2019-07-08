function [OsciFrac] = IRASAPower_EEG_wrapper(spec,oscifrac,frange,normalize,method)

OsciFrac = zeros(1,size(spec.osci,2));

if nargin < 5
    method = 'mean';
end

if nargin < 4
    normalize = 'no';
end

if nargin < 3
    frange = [0.5 50];
end

disp(' ')
if strcmpi(oscifrac,'osci')
    disp('Computing oscillatory power...')
else
    disp('Computing fractal power...')
end

freqindex = intersect(find(spec.freq > frange(1)),find(spec.freq < frange(2)));
freqs = spec.freq(freqindex);

for c = 1:size(spec.frac,2)
    fprintf([num2str(c) ' ']);
    %osci = sgolayfilt(spec.osci(:,c),5,1501);
    %osci = osci(freqindex);
    %frac = sgolayfilt(spec.frac(:,c),5,1501);
    pow = spec.(oscifrac)(freqindex,c);
    
    if strcmpi(method,'trapz')
        OsciFrac(c) = trapz(freqs,pow);
    elseif strcmpi(method,'mean')
        OsciFrac(c) = nanmean(pow);
    end
    if strcmpi(normalize,'mixd')
        if strcmpi(method,'trapz')
            OsciFrac(c) = OsciFrac(c)/trapz(freqs,spec.mixd(freqindex,c));
        elseif strcmpi(method,'mean')
            OsciFrac(c) = OsciFrac(c)/nanmean(spec.mixd(freqindex,c));
        end
    elseif iscell(normalize)
        newfreqs = intersect(find(spec.freq > normalize{2}(1)),find(spec.freq < normalize{2}(2)));
        if strcmpi(method,'trapz')
            OsciFrac(c) = OsciFrac(c)/trapz(spec.freq(newfreqs),spec.(normalize{1})(newfreqs,c));
        elseif strcmpi(method,'mean')
            OsciFrac(c) = OsciFrac(c)/nanmean(spec.(normalize{1})(newfreqs,c));
        end
    end
end