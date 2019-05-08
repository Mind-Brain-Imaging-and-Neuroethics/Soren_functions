function [DFA_irasa] = DFA_IRASA_EEG_handle(spec)

DFA_irasa = zeros(1,size(spec.osci,2));

disp(' ')
disp('Computing DFA from IRASA spectrum...')

amp = sqrt(2*spec.frac);

phases = 2*pi*rand(length(amp),1);

freqSig = amp.*e^(j*phases);

timeSeries = ifft(freqSig);

freqindex = intersect(find(spec.freq > 0.1),find(spec.freq < 70));
freqs = spec.freq(freqindex);
    
for c = 1:size(spec.osci,2)
    fprintf([num2str(c) ' ']);
    osci = sgolayfilt(spec.osci(:,c),5,1501);
    osci = osci(freqindex);
    frac = sgolayfilt(spec.frac(:,c),5,1501);
    frac = frac(freqindex);
    
    DFA_irasa(c) = trapz(freqs,osci)/trapz(freqs,frac);
end