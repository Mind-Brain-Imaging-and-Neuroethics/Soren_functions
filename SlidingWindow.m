function [fwindows,tpoints,meanwindow] = SlidingWindow(funchandle,data,winsize,overlap)

ii=1;   % Windows counter
while true
    % Begining and ending of the current window
    SWindow=[round(1+(ii-1)*winsize*(1-overlap)), round(winsize+(ii-1)*winsize*(1-overlap))];
    % Check if index exceeds vector dimensions. If so, break!
    if SWindow(2)>=length(data)
        break
    end
    % ACF computation into the window (normalized between -1 and 1)
    fwindows(ii) = funchandle(data(SWindow(1):SWindow(2)));
    tpoints(ii) = SWindow(1);
    % Next window
    ii=ii+1;
end
meanwindow = mean(fwindows);