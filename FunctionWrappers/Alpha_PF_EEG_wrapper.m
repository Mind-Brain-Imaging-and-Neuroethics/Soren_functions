function [pfout] = Alpha_PF_EEG_wrapper(EEG)

disp('Computing alpha peak frequency...')
disp('')

[psum] = restingIAF(EEG.data,EEG.nbchan,3,[1 40],EEG.srate,[5 15],11,5);

if ~isnan(psum.paf)
    pf = psum.paf;
else
    pf = psum.cog;
end

pfout = repmat(pf,1,EEG.nbchan);
