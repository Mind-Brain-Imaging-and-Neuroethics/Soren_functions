function [indices,allindices,flatindices] = Ivar_findConds(conditions,conditionnames)

indices = cell(1,2);
for c = 1:length(conditions)
   if iscell(conditions{c})
       indices{c} = cell(1,length(conditions{c}));
   end
end

allindices = [];
for q = 1:length(conditions) %for each condition
    for c = 1:length(conditionnames) %loop through each file
        skipFile = 0;
        if ~iscell(conditions{q})
            for cc = 1:length(conditions{q}) %loop through each character in the conditions
                %if ~(conditionnames{c}{1}(cc) == conditions{q}(cc)) && ~(conditions{q}(cc) == 'x')
                if ~(conditionnames{c}(cc) == conditions{q}(cc)) && ~(conditions{q}(cc) == 'x')
                    skipFile = 1;
                end
            end
            if ~skipFile
                indices{q} = [indices{q} c];
                allindices = [allindices c];
            end
            
        else
            for i = 1:length(conditions{q}) %for combining conditions
                skipFile = 0;
                for cc = 1:length(conditions{q}{i}) %loop through each character in the conditions
                    %if ~(conditionnames{c}{1}(cc) == conditions{q}{i}(cc)) && ~(conditions{q}{i}(cc) == 'x')
                    if ~(conditionnames{c}(cc) == conditions{q}{i}(cc)) && ~(conditions{q}{i}(cc) == 'x')
                        skipFile = 1;
                    end
                end
                if ~skipFile
                    indices{q}{i} = [indices{q}{i} c];
                    allindices = [allindices c];
                end
            end
        end
            
    end
end

flatindices = cell(1,length(conditions));
if nargout > 2
    for c = 1:length(flatindices)
        if iscell(indices{c})
            for cc = 1:length(indices{c})
                flatindices{c} = horzcat(flatindices{c},indices{c}{cc});
            end
        else
            flatindices{c} = indices{c};
        end
    end
end


