function TTV_ERSP_FDR_correct(settings)

load([settings.outputdir '/' settings.datasetname '_results.mat'])
if isfield(settings,'rest')
    load([settings.outputdir '/' settings.datasetname '_restmeas.mat'])
end

fields = alloutputs.fdrfields;

allpvals = [];
sourceindex = [];
for c = 1:length(fields)
    data = getfield_nest(alloutputs,fields{c});
    for cc = 1:length(data)
        if isfield(data{cc},'posclusters') && ~isempty(data{cc}.posclusters)
            allpvals = [allpvals horz(extractfield(data{cc}.posclusters,'prob'))];
            for ccc = 1:length(extractfield(data{cc}.posclusters,'prob'))
                sourceindex = [sourceindex {['alloutputs.' fields{c} '{' num2str(cc) '}.posclusters(' num2str(ccc) ').prob']}];
            end
        end
        if isfield(data{cc},'negclusters') && ~isempty(data{cc}.negclusters)
            allpvals = [allpvals horz(extractfield(data{cc}.negclusters,'prob'))];
            for ccc = 1:length(extractfield(data{cc}.negclusters,'prob'))
                sourceindex = [sourceindex {['alloutputs.' fields{c} '{' num2str(cc) '}.negclusters(' num2str(ccc) ').prob']}];
            end
        end
    end
end

if isfield(settings,'rest')
    fields = restmeas.fdrfields;
    
    for c = 1:length(fields)
        data = getfield_nest(restmeas,fields{c});
        for cc = 1:length(data)
            if isfield(data{cc},'posclusters') && ~isempty(data{cc}.posclusters)
                allpvals = [allpvals horz(extractfield(data{cc}.posclusters,'prob'))];
                for ccc = 1:length(extractfield(data{cc}.posclusters,'prob'))
                    sourceindex = [sourceindex {['restmeas.' fields{c} '{' num2str(cc) '}.posclusters(' num2str(ccc) ').prob']}];
                end
            end
            if isfield(data{cc},'negclusters') && ~isempty(data{cc}.negclusters)
                allpvals = [allpvals horz(extractfield(data{cc}.negclusters,'prob'))];
                for ccc = 1:length(extractfield(data{cc}.negclusters,'prob'))
                    sourceindex = [sourceindex {['restmeas.' fields{c} '{' num2str(cc) '}.negclusters(' num2str(ccc) ').prob']}];
                end
            end
            if isfield(data{cc},'sobel')
                allpvals = [allpvals data{cc}.sobel.p];
                sourceindex = [sourceindex {['restmeas.' fields{c} '{' num2str(cc) '}.sobel.p']}];
                allpvals = [allpvals data{cc}.montecarlo.p];
                sourceindex = [sourceindex {['restmeas.' fields{c} '{' num2str(cc) '}.montecarlo.p']}];
            end
        end
    end
end

newpvals = mafdr(allpvals,'BHFDR',true);

for c = 1:length(sourceindex)
    p = newpvals(c);
    eval([sourceindex{c} ' =  p;']);
end

save([settings.outputdir '/' settings.datasetname '_results.mat'],'alloutputs')
if isfield(settings,'rest')
    save([settings.outputdir '/' settings.datasetname '_restmeas.mat'],'restmeas')
end

