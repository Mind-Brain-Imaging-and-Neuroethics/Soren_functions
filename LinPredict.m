function [predWeights,predValue] = LinPredict(timeSeries,rx,PLE)
%currently doesn't support predWeights for fbm sequences

K = length(timeSeries);

if length(timeSeries) == 1
    predWeights = 1;
    predValue = timeSeries(1);
else
    if PLE < 1
        %this only works for a zero-meaned time series, so subtract the mean
        origmean = mean(timeSeries);
        zmtimeSeries = timeSeries-origmean;
        
        %V = var(zmtimeSeries);

        %rx = CovarFromHurst(PLE,zmtimeSeries,V);
        
        Rk = zeros(K);
        
        for p = 1:K
            for pp = 1:K
                %Rk(p,pp) = CovarFromHurst(PLE,zmtimeSeries,V,abs(p-pp));
                Rk(p,pp) = rx(abs(p-pp)+1);
            end
        end
        %Rk = Rk + Rk';
        %Rk = Rk./(eye(K)+1);

        
        if isinf(cond(Rk))
            predWeights = zeros(1,length(timeSeries));
            predValue = mean(timeSeries);
        else
            Rk = inv(Rk);
            predWeights = Rk*rx(1:K)';
            
            predValue = sum(predWeights'.*fliplr(zmtimeSeries))+origmean;
        end
        
    else
        
        fgn = zeros (1,K-1);
        for c = 1:K-1
            fgn(c) = timeSeries(c+1) - timeSeries(c);
        end
        zmfgn = fgn - mean(fgn);
        
        rx = CovarFromHurst(PLE-2,zmfgn);
        
        K = length(zmfgn);
        
        Rk = zeros(K);
        
        for p = 1:K
            for pp = 1:K
                Rk(p,pp) = CovarFromHurst(PLE-2,zmfgn,abs(p-pp));
            end
        end
        
        
        if isinf(cond(Rk))
            predValueFGN = mean(fgn);
        else
            Rk = inv(Rk);
            
            predWeightsFGN = Rk*rx';
            
            predValueFGN = sum(predWeightsFGN.*zmfgn')+mean(fgn);
        end
        
        predValue = timeSeries(end) + predValueFGN;
        
        predWeights = zeros(1,length(timeSeries));
        
        %predWeights =
    end
end

binWidth = 1/16;
predValue = round(predValue/binWidth)*binWidth;

    
