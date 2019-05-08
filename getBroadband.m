function [dataout] = getBroadband(datain)
    dataout = [];
    for c = 1:size(datain,1)
        if mod(c,5) == 1
        dataout = vertcat(dataout,datain(c,:));
        end
    end
end