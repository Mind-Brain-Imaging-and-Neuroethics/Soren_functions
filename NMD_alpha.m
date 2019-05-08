function [dataout] = NMD_alpha(datain,srate,bandpass)

parfor c = 1:size(datain,1)
    currdata = diff([0 datain(c,:)]);
    for cc = 1:10
        diffmode = nmd_mod(currdata,srate,'S',0,'S0',0,'TFRtype','WT','fmin',bandpass(1),'fmax',bandpass(2),'ModeNum',1);
        tmp2 = diff(angle(hilbert(diffmode)))*srate/(2*pi);
        if median(tmp2) > 8 && median(tmp2) < 14
            tmp = cumsum(diffmode);
            mode = nmd_mod(tmp,srate,'S',0,'S0',0,'TFRtype','WT','fmin',bandpass(1),'fmax',bandpass(2),'ModeNum',1);
            dataout(c,:) = mode;
            break
        else
            currdata = currdata - diffmode;
        end
    end
end