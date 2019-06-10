function [ROI] = Plot_parc_data(datain,sourcemodel,atlas,varargin)

%if size(corrdata)

    if CheckInput(varargin,'mask')
        roimask = EasyParse(varargin,'mask');
    end

for c = 1:length(sourcemodel.pos)
    if isfield(atlas,'parcellation')
        plotdata(c) = datain(atlas.parcellation(c));
    elseif isfield(atlas,'parcels')
        plotdata(c) = datain(atlas.parcels(c));
    end
    if CheckInput(varargin,'mask')
        if isfield(atlas,'parcellation')
       plotmask(c) = roimask(atlas.parcellation(c)); 
        else
            plotmask(c) = roimask(atlas.parcels(c));
        end
    end
end


if ~EasyParse(varargin,'SinglePlot','on')
    subplot(2,2,[1 3])
    bnd.pnt = sourcemodel.pos;
    bnd.tri = sourcemodel.tri;
    
    ft_plot_mesh(bnd,'vertexcolor',plotdata');
    clim = get(gca,'CLim');
    colorbar
    
    cort_size = length(find(sourcemodel.brainstructure == 1));
    
    subplot(2,2,2)
    bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 1),:);
    trindx = find(max(sourcemodel.tri,[],2) <= cort_size);
    bnd.tri = sourcemodel.tri(trindx,:);
    
    set(gca,'CLim',clim);
    
    subplot(2,2,4)
    bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 2),:);
    trindx = find(min(sourcemodel.tri,[],2) >= cort_size);
    bnd.tri = sourcemodel.tri(trindx,:);
    bnd.tri = bnd.tri - cort_size;
    
    ft_plot_mesh(bnd,'vertexcolor',plotdata((length(plotdata)/2+1):end)','edgealpha',0);
    set(gca,'CLim',clim)
else
    bnd.pnt = sourcemodel.pos;
    bnd.tri = sourcemodel.tri;
    if CheckInput(varargin,'mask')
        plotmask(find(~plotmask)) = NaN;
        plotdata_masked = plotdata.*plotmask;
        ft_plot_mesh(bnd,'vertexcolor',plotdata_masked');
        hold on
        ft_plot_mesh(bnd,'vertexcolor',plotdata','facealpha',0.5)
    else
        ft_plot_mesh(bnd,'vertexcolor',plotdata');
    end
end

