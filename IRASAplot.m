function IRASAplot(spec,bandpass)

if nargin < 2
   bandpass = [spec.freq(1)-1 spec.freq(end)+1]; 
end
freqs = intersect(find(spec.freq > bandpass(1)), find(spec.freq < bandpass(2)));
% if strcmpi(oscifrac,'osci')
%     plot(spec.freq(freqs),spec.(oscifrac)(freqs));
% else
%     loglog(spec.freq(freqs),spec.(oscifrac)(freqs));
% end

subplot(2,2,[1 2])
loglog(spec.freq(freqs),mean(spec.mixd(freqs,:),2))
subplot(2,2,3)
plot(spec.freq(freqs),mean(spec.osci(freqs,:),2))
subplot(2,2,4)
loglog(spec.freq(freqs),mean(spec.frac(freqs,:),2))
