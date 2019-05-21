function [pp1,power_data,power_freq] = JF_power_law(time_series,TR,low_range,high_range,varargin)
%     HS = spectrum.welch;
Fs = 1/TR;
%     nfft = 2^nextpow2(length(time_series));

[pxx,f] = pwelch(time_series,[],[],2^nextpow2((3/low_range)*Fs),Fs); %want 3 cycles of lowest frequency in window
%     power_spec = psd(HS,time_series,'NFFT',nfft,'Fs',Fs);
power_data = pxx;
power_freq = f;

slope_index = find(power_freq > low_range & power_freq < high_range);
power_freq = linspace(min(f(slope_index)),max(f(slope_index)),length(f(slope_index)));
p = polyfit(log(power_freq)',log(power_data(slope_index)),1);
pp1 = -p(1);


if length(varargin) > 0 && any(varargin{1} == 'Plot')
    y = p(2) + p(1)*log(f(slope_index));
    loglog(f(slope_index),pxx(slope_index));
    hold on;
%    loglog(f(slope_index),exp(y),'r--');
    xlabel('Log Frequency')
    ylabel('Log Power')
    title(['Estimated PLE is ' num2str(pp1)])
end


%[b,bint,r,rint,stats] = regress(log(power_data(slope_index)),[ones(length(power_freq(slope_index)),1) log(power_freq')*p(1)+ p(2)]);
%r_resi = stats(4);
