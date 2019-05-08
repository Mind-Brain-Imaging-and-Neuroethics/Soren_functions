function [data] = ft_autoartifact(data,types,channels,interactive)

if nargin < 2
    interactive = 0;
end

if ~isempty(find(strcmpi(types,'jump')))
    % Jump artifacts
    cfg = [];
    
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel = channels;
    cfg.artfctdef.zvalue.cutoff = 20;
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 0;
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.cumulative = 'yes';
    cfg.artfctdef.zvalue.medianfilter = 'yes';
    cfg.artfctdef.zvalue.medianfiltord = 9;
    cfg.artfctdef.zvalue.absdiff = 'yes';
    
    % % make the process interactive
    if interactive
        cfg.artfctdef.zvalue.interactive = 'yes';
    end
    
    [cfg, artifact_jump] = ft_artifact_zvalue(cfg,data);
end

if ~isempty(find(strcmpi(types,'muscle')))
    % Muscle artifacts
    cfg            = [];
    
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel      = channels;
    cfg.artfctdef.zvalue.cutoff       = 4;
    cfg.artfctdef.zvalue.trlpadding   = 0;
    cfg.artfctdef.zvalue.fltpadding   = 0;
    cfg.artfctdef.zvalue.artpadding   = 0.1;
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.bpfilter     = 'yes';
    cfg.artfctdef.zvalue.bpfreq       = [110 140];
    cfg.artfctdef.zvalue.bpfiltord    = 9;
    cfg.artfctdef.zvalue.bpfilttype   = 'but';
    cfg.artfctdef.zvalue.hilbert      = 'yes';
    cfg.artfctdef.zvalue.boxcar       = 0.2;
    
    % make the process interactive
    if interactive
        cfg.artfctdef.zvalue.interactive = 'yes';
    end
    
    [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,data);
end

if ~isempty(find(strcmpi(types,'EOG')))
    % EOG
    cfg            = [];
    
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel     = 'EOG';
    cfg.artfctdef.zvalue.cutoff      = 4;
    cfg.artfctdef.zvalue.trlpadding  = 0;
    cfg.artfctdef.zvalue.artpadding  = 0.1;
    cfg.artfctdef.zvalue.fltpadding  = 0;
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.bpfilter   = 'yes';
    cfg.artfctdef.zvalue.bpfilttype = 'but';
    cfg.artfctdef.zvalue.bpfreq     = [2 15];
    cfg.artfctdef.zvalue.bpfiltord  = 4;
    cfg.artfctdef.zvalue.hilbert    = 'yes';
    
    % feedback
    if interactive
        cfg.artfctdef.zvalue.interactive = 'yes';
    end
    
    [cfg, artifact_EOG] = ft_artifact_zvalue(cfg,data);
end


cfg=[];
cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
if exist('artifact_EOG','var')
    cfg.artfctdef.eog.artifact = artifact_EOG; 
end
if exist('artifact_jump','var')
    cfg.artfctdef.jump.artifact = artifact_jump;
end
if exist('artifact_muscle','var')
    cfg.artfctdef.muscle.artifact = artifact_muscle;
end


data = ft_rejectartifact(cfg,data);