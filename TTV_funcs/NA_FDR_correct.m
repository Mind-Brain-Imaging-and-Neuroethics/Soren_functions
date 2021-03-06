function NA_FDR_correct(settings)

load([settings.outputdir '/' settings.datasetname '_results.mat'])
if isfield(settings,'rest')
    load([settings.outputdir '/' settings.datasetname '_restmeas.mat'])
end

fields = alloutputs.fdrfields;

allpvals = [];
sourceindex = [];
for c = 1:length(fields)
    data = getfield_nest(alloutputs,fields{c});
    if iscell(data)
        data = reshape(data,[],1);
        for cc = 1:length(data)
            if isfield(data{cc},'posclusters') && ~isempty(data{cc}.posclusters)
                tmp = horz(extractfield(data{cc}.posclusters,'prob'));
                allpvals = [allpvals tmp(1)];
                %for ccc = 1:length(extractfield(data{cc}.posclusters,'prob'))
                ccc = 1;
                    sourceindex = [sourceindex {['alloutputs.' fields{c} '{' num2str(cc) '}.posclusters(' num2str(ccc) ').prob']}];
                %end
            end
            if isfield(data{cc},'negclusters') && ~isempty(data{cc}.negclusters)
                tmp = horz(extractfield(data{cc}.negclusters,'prob'));
                allpvals = [allpvals tmp(1)];
                %for ccc = 1:length(extractfield(data{cc}.negclusters,'prob'))
                ccc = 1;    
                sourceindex = [sourceindex {['alloutputs.' fields{c} '{' num2str(cc) '}.negclusters(' num2str(ccc) ').prob']}];
                %end
            end
            if ~(isfield(data{cc},'posclusters') && ~isempty(data{cc}.posclusters)) && ...
                ~(isfield(data{cc},'negclusters') && ~isempty(data{cc}.negclusters))
                allpvals = [allpvals 1];
            end
        end
    elseif isempty(data)
        %do nothing
    elseif isnumeric(data)
        allpvals = [allpvals horz(data)];
        sourceindex = [sourceindex {['alloutputs.' fields{c}]}];
    end
end

if isfield(settings,'rest')
    fields = restmeas.fdrfields;
    
    for c = 1:length(fields)
        data = getfield_nest(restmeas,fields{c});
        if iscell(data)
            data = reshape(data,[],1);
            for cc = 1:length(data)
                if isfield(data{cc},'posclusters') && ~isempty(data{cc}.posclusters)
                    tmp = horz(extractfield(data{cc}.posclusters,'prob'));
                    allpvals = [allpvals tmp(1)];
                    %for ccc = 1:length(extractfield(data{cc}.posclusters,'prob'))
                    ccc = 1;
                    sourceindex = [sourceindex {['restmeas.' fields{c} '{' num2str(cc) '}.posclusters(' num2str(ccc) ').prob']}];
                    %end
                end
                if isfield(data{cc},'negclusters') && ~isempty(data{cc}.negclusters)
                    tmp = horz(extractfield(data{cc}.negclusters,'prob'));
                    allpvals = [allpvals tmp(1)];
                    %for ccc = 1:length(extractfield(data{cc}.negclusters,'prob'))
                    ccc = 1;
                    sourceindex = [sourceindex {['restmeas.' fields{c} '{' num2str(cc) '}.negclusters(' num2str(ccc) ').prob']}];
                    %end
                end
                if isfield(data{cc},'sobel')
                    allpvals = [allpvals data{cc}.sobel.p];
                    sourceindex = [sourceindex {['restmeas.' fields{c} '{' num2str(cc) '}.sobel.p']}];
                    allpvals = [allpvals data{cc}.montecarlo.p];
                    sourceindex = [sourceindex {['restmeas.' fields{c} '{' num2str(cc) '}.montecarlo.p']}];
                end
                if ~(isfield(data{cc},'posclusters') && ~isempty(data{cc}.posclusters)) && ...
                    ~(isfield(data{cc},'negclusters') && ~isempty(data{cc}.negclusters)) && ...
                    ~isfield(data{cc},'sobel')
                    allpvals = [allpvals 1];
                end
            end
        elseif isnumeric(data)
            allpvals = [allpvals horz(data)];
            sourceindex = [sourceindex {['alloutputs.' fields{c}]}];
        end
        
    end
end

newpvals = mafdr(allpvals,'BHFDR',true);

for c = 1:length(sourceindex)
    p = newpvals(c);
    eval([sourceindex{c} ' =  p;']);
end

save([settings.outputdir '/' settings.datasetname '_results_FDR.mat'],'alloutputs')
if isfield(settings,'rest')
    save([settings.outputdir '/' settings.datasetname '_restmeas_FDR.mat'],'restmeas')
end

