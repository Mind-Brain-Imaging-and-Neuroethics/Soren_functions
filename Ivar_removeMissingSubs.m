function [removesubs] = Ivar_removeMissingSubs(indices,subnames,ALLEEG)

deletedsubs = {};
    removesubs = [];
    
    for q = 1:2
        if iscell(indices{q})
            notq = mod(q,2)+1;
            for c = 1:length(indices{notq})
                sub1 = subnames{indices{q}{1}(c)};
                findsub = find(strcmpi(sub1,subnames(indices{q}{2})));
                set1 = find(allindices == indices{q}{1}(c));
                
                if ~isempty(findsub)
                    set2 = find(allindices == indices{2}{2}(findsub));
                    ALLEEG(set1).data = cat(3,ALLEEG(set1).data,ALLEEG(set2).data);
                    ALLEEG(set1) = eeg_checkset(ALLEEG(set1));
                    removesubs = [removesubs set2];
                else
                    removesubs = [removesubs set1];
                    deletedsubs = [deletedsubs {sub1}];
                end
            end
        end
        
        if ~isempty(deletedsubs)
            for c = 1:length(deletedsubs)
                if ~iscell(indices{notq})
                    for cc = 1:length(indices{notq})
                        if strcmpi(deletedsubs{c},subnames{indices{notq}(cc)})
                            removesubs = [removesubs find(allindices == indices{1}(cc))]
                        end
                    end
                else
                    for ccc = 1:2
                        for cc = 1:length(indices{notq}{ccc})
                            if strcmpi(deletedsubs{c},subnames{indices{notq}{ccc}(cc)})
                                removesubs = [removesubs find(allindices == indices{1}{ccc}(cc))]
                            end
                        end
                    end
                end
            end
        end
    end