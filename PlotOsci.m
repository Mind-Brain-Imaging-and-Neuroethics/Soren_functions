function PlotOsci(spec,frange,electrode)

figure

freqs = intersect(find(spec.freq > frange(1)),find(spec.freq < frange(2)));

osciplot = sgolayfilt(spec.osci(:,electrode),11,1501);

plot(spec.freq(freqs),osciplot(freqs))