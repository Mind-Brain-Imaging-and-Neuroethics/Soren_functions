function IRASAplot(spec,oscifrac,bandpass)

freqs = intersect(find(spec.freq > bandpass(1)), find(spec.freq < bandpass(2)));
if strcmpi(oscifrac,'osci')
    plot(spec.freq(freqs),spec.(oscifrac)(freqs));
else
    loglog(spec.freq(freqs),spec.(oscifrac)(freqs));
end