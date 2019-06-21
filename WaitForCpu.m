function WaitForCpu(targUsage,command,varargin)

% if cpu usage is consistently below the target, run the command

varargin = setdefault(varargin,'historylength',6);
varargin = setdefault(varargin,'checkinterval',300);

avgusage = ones(1,EasyParse(varargin,'historylength'))*100;
checkinterval = EasyParse(varargin,'checkinterval');

while 1
    [~,usage] = system('top -bn2 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk ''{print 100 - $1}''');
    usage = str2num(usage);
    fprintf([num2str(usage(2)) newline])
    avgusage = avgusage([2:end 1]);
    avgusage(end) = usage(2);
    
    if isempty(find(avgusage > targUsage))
       eval(command);
       break;
    end
    
    pause(checkinterval) % pause before checking again
end

end
