function [output] = MakeFGN(PLE,seqLen)
%takes a PLE, NOT a hurst exponent as an argument

Hexpfbm = (PLE+1)/2; %converts the input PLE to a hurst exponent for a corresponding fractional brownian motion sequence

if seqLen < 1000
fbmseq = wfbm(Hexpfbm,1000);
else
    fbmseq = wfbm(Hexpfbm,seqLen+1);
end

output = diff(fbmseq);



