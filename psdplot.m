function psdplot(data,fs,range)
% Quick function for plotting power spectra

if size(data,2) == 1
    data = data';
end

if nargin > 2
    for c = 1:size(data,1)
        [pxx(c,:),f(c,:)] = pwelch(data(c,:),[],[],2^nextpow2((3/range(1))*fs),fs); %want 3 cycles of lowest frequency in window
        frange = intersect(find(f(c,:) > range(1)),find(f(c,:) < range(2)));
        loglog(f(c,frange),pxx(c,frange))
        hold on
    end
else
    for c = 1:size(data,1)
        [pxx(c,:),f(c,:)] = pwelch(data(c,:),[],[],[],fs);
        loglog(f(c,2:end),pxx(c,2:end))
        hold on
    end
end




