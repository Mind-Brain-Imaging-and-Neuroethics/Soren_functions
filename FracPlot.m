function FracPlot(spec_orig,spec_new,frange,chanindx)

if nargin < 4
   chanindx = [1 1];
end

findex_orig = intersect(find(spec_orig.freq > frange(1)),find(spec_orig.freq < frange(2)));
findex_new = intersect(find(spec_new.freq > frange(1)),find(spec_new.freq < frange(2)));

if ischar(chanindx) && strcmpi(chanindx,'mean')
   spec_orig.osci = mean(spec_orig.osci,2);
   spec_orig.frac = mean(spec_orig.frac,2);
   spec_new.osci = mean(spec_new.osci,2);
   spec_new.frac = mean(spec_new.frac,2);
   chanindx = [1 1];
end
clf
subplot(2,1,1)
%loglog(spec_orig.freq(findex_orig),spec_orig.osci(findex_orig,chanindx(1)),'b')
plot(spec_orig.freq(findex_orig),spec_orig.osci(findex_orig,chanindx(1)),'b')
hold on
%loglog(spec_new.freq(findex_new),spec_new.osci(findex_new,chanindx(2)),'r')
plot(spec_new.freq(findex_new),spec_new.osci(findex_new,chanindx(2)),'r')
legend({'Original spectrum','Spectrum - mode(s)'})
title('Oscillatory component')
xlabel('Frequency (Hz)')
ylabel('Power')

subplot(2,1,2)
loglog(spec_orig.freq(findex_orig),spec_orig.frac(findex_orig,chanindx(1)),'b')
%plot(spec_orig.freq(findex_orig),spec_orig.frac(findex_orig,chanindx(1)),'b')
hold on
loglog(spec_new.freq(findex_new),spec_new.frac(findex_new,chanindx(2)),'r')
%plot(spec_new.freq(findex_new),spec_new.frac(findex_new,chanindx(2)),'r')
legend({'Original spectrum','Spectrum - mode(s)'})
title('Fractal component')
xlabel('Frequency (Hz)')
ylabel('Power')