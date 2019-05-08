function [dataout] = Mean_downsample(datain,N)


for c = 1:ceil(length(datain)/N)
    if c*N < length(datain)
        dataout(c) = mean(datain(((c-1)*N+1):(c*N)));
    else
        dataout(c) = mean(datain(((c-1)*N+1):end));
    end
end