function WaitForCpu(targUsage,command)

% keep a 18-minute moving average of usage - if this is consistently below
% the target, run the command

avgusage = ones(1,6)*100;
while 1
    usage = system('top -b -d1 -n1|grep -i "Cpu(s)"|head -c21|cut -d '' '' -f2|cut -d ''%'' -f1')
    avgusage = avgusage([2:end 1]);
    avgusage(end) = usage;
    
    if isempty(find(avgusage > targUsage,1))
       eval(command);
       break;
    end
    
    pause(3) % pause for 3 mins before checking again
end

end