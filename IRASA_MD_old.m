function [modes,fractal,specs,unusedmode,unusedspec] = IRASA_MD(datain,srate,bandpass,varargin)

count = 1;
specs(count) = amri_sig_fractal(datain,srate);
criterion = 0;

while ~criterion
    if count == 1
       newdata = datain; 
    end
    olddata = newdata;
    [modes(count,:)] = nmd(double(olddata),srate,'ModeNum',1,'S',2,'Template','Fast');
    newdata = newdata - modes(count,:);
    specs(count+1) = amri_sig_fractal(newdata,srate);
    realfreqs = intersect(find(specs(count+1).freq > bandpass(1)),find(specs(count+1).freq < bandpass(2)));
    norm_osci = norm(specs(count).osci(realfreqs) - specs(count+1).osci(realfreqs));
    norm_frac = norm(specs(count).frac(realfreqs) - specs(count+1).frac(realfreqs));
    if norm_osci < 2*norm_frac
       unusedmodes = modes(count,:);
       modes(count,:) = [];
       unusedspec = specs(count+1);
       specs(count+1) = [];
       fractal = olddata;
       criterion = 1;
    else
        count = count+1;
    end
end

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