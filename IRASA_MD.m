function [modes,oscireduce,fracreduce] = IRASA_MD(datain,srate,bandpass,varargin)

count = 1;
specs(1) = amri_sig_fractal(datain,srate);
freqs = intersect(find(specs(1).freq > bandpass(1)),find(specs(1).freq < bandpass(2)));
criterion = 0;

%[modes] = nmd(double(datain),srate,'ModeNum',5,'S',2,'Template','fast');

for c = 1:5
   tmpdata = datain - modes(c,:);
   specs(c+1) = amri_sig_fractal(tmpdata,srate);
   oscireduce(c) = norm(specs(c+1).osci(freqs) - specs(1).osci(freqs));
   %fraccorr(c) = norm(specs(c+1).frac(freqs) - specs(1).frac(freqs));
    fraccorr(c) = corr(specs(c+1).frac(freqs),specs(1).frac(freqs));
   %    figure 
%    subplot(2,1,1)
%    plot(specs(1).freq(freqs),specs(1).osci(freqs))
%    hold on
%    plot(specs(1).freq(freqs),specs(c+1).osci(freqs))
%     subplot(2,1,2)
%        loglog(specs(1).freq(freqs),specs(1).frac(freqs))
%    hold on
%    loglog(specs(1).freq(freqs),specs(c+1).frac(freqs))
end







% while ~criterion
%     if count == 1
%        newdata = datain; 
%     end
%     olddata = newdata;
%     newdata = newdata - modes(count,:);
%     specs(count+1) = amri_sig_fractal(newdata,srate);
%     realfreqs = intersect(find(specs(count+1).freq > bandpass(1)),find(specs(count+1).freq < bandpass(2)));
%     norm_osci = norm(specs(count).osci(realfreqs) - specs(count+1).osci(realfreqs));
%     norm_frac = norm(specs(count).frac(realfreqs) - specs(count+1).frac(realfreqs));
%     if norm_osci < 2*norm_frac
%        unusedmodes = modes(count,:);
%        modes(count,:) = [];
%        unusedspec = specs(count+1);
%        specs(count+1) = [];
%        fractal = olddata;
%        criterion = 1;
%     else
%         count = count+1;
%     end
% end

% if EasyParse(varargin,'Plot','on')
%    n = size(modes,1);
%    f = factor(n+1);
%    if length(f) < 2
%       f = [1 f]; 
%    end
%    
%    figure
%    for c = 1:(n-1)
%       subplot(f(1),f(2),c)
%       plot(specs(c).freq(realfreqs),specs(c).osci(realfreqs));
%    end
%    subplot(f(1),f(2),n)
%    for c = 1:(n-1)
%       plot(specs(c).freq(realfreqs),specs(c).osci(realfreqs));
%       hold on
%    end
%    
%    
%    figure
%       for c = 1:(n-1)
%       subplot(f(1),f(2),c)
%       plot(specs(c).freq(realfreqs),specs(c).frac(realfreqs));
%    end
%    subplot(f(1),f(2),n)
%    for c = 1:(n-1)
%       plot(specs(c).freq(realfreqs),specs(c).frac(realfreqs));
%       hold on
%    end
%    
% end