function [PLE] = PedroPLE(data)



%%%%Loop through all the channels within each file and do a windowed FFT
for c = 1:size(data,1)
    tmpdata = data(c,:);
    
    %Determine the size of our data for further processing
    nSamples = length(tmpdata);
    
    %Define the window length in which we want to divide the continuous file
    winSize = 2000;
    
    %Overlap between windows, in percentage
    OverlapStep = 50;
    
    %Use overlap size to pre-arrange the data
    if OverlapStep > 0
        Overlap = floor((OverlapStep*winSize) / 100);
        nFrames=floor(nSamples/Overlap)-1;
    else
        Overlap= winSize;
        nFrames=floor(nSamples/Overlap)-1;
    end
    
    %Go through all the created windows and extract one FFT per each before averaging them
    k = 1;
    z = 1;
    while ( (k+winSize-1) <= nSamples )
        z;
        FrameSignal = tmpdata(k:k+winSize-1);
        
        w = flattopwin(winSize, 'periodic'); %Type of windowing to do (Hanning, Flattop, etc...we choose FlatTop because empirically it distorts less the estimations of magnitude (amplitude))
        wdata = FrameSignal(:).*w; %Convolute our epochs with the window we want
        NFFT = 2^nextpow2(length(wdata)); %Calculate size of the FFT as a power of 2 to improve performance
        
        %Extract our FFTs and average on the fly the amplitude
        if z == 1
            Y = abs(fft(FrameSignal, NFFT)/2);
        else
            Y = (Y+abs(fft(FrameSignal, NFFT)/2))/2;
        end
        k=k+Overlap;
        z = z+1;
    end
    
    %Save the FFTs amplitudes for each subject for each channel
    fftdata(c,:) = Y;
    
    %Create the vector of our frequencies, based on the sampling rate (500Hz) and the size of our FFT, so we can know where we are - adjust accordingly
    f = 500*(0:(NFFT/2))/NFFT;
    
    %use MATLAB's fit() function to do a linear fit between the log transformation of both the frequency range and the spectrum
    p = polyfit(log(f(2:185))', log(Y((2:185)))', 1);
    
    %the fit is done using a y=mx + b equation (1st degree polynomial, i.e., linear fit), where m is the slope - the variable fitobject contains the various parameters of the fit
    PLE(c) = abs(p(1,1));
    
end
