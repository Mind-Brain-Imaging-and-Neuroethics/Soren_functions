function [meanPLEirasa,PLEsirasa,specobject] = IRASA_window(EEGdata,srate,overlap)

winsize = 120;
numwins = floor(length(EEGdata)/(srate*winsize*overlap))-1;


PLEsraw = zeros(size(EEGdata,1),numwins);

for c = 1:numwins
    tmp(c) = amri_sig_fractal(EEGdata(:,((c-1)*(srate*winsize*overlap)+1):(c-1)*(srate*winsize*overlap)+srate*winsize),srate);
    a = amri_sig_plawfit(tmp(c),[0.1 70]);
    PLEsirasa(:,c) = a.Beta;
end

for c = 1:size(EEGdata,1)
   meanPLEirasa(c) = nan_mean(PLEsirasa(c,:)); 
end

if nargout > 1
    specobject = struct;
    
    tmp2 = struct;
    tmp2.Freq = [];
    tmp2.Frac = [];
    tmp2.Osci = [];
    
    for c = 1:2
        tmp2.Freq = cat(3,tmp2.Freq,tmp(c).freq);
        tmp2.Frac = cat(3,tmp2.Frac,tmp(c).frac);
        tmp2.Osci = cat(3,tmp2.Osci,tmp(c).osci);
    end
    
    specobject.Freq = nan_mean(tmp2.Freq,3);
    specobject.Frac = nan_mean(tmp2.Frac,3);
    specobject.Osci = nan_mean(tmp2.Osci,3);
end



