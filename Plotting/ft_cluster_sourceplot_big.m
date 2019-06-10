function [ROI] = ft_cluster_sourceplot(datain,sourcemodel,atlas,roimask,varargin)


if ~CheckInput(varargin,'colormaps')
    cmap1= parula;
    cmap2 = gray;
else
    cmaps = EasyParse(varargin,'colormaps')
    cmap1 = cmaps{1};
    cmap2 = cmaps{2};
end
cindex = linspace(min(datain),max(datain),64);

for c = 1:length(sourcemodel.pos)
    if isfield(atlas,'parcellation')
        plotdata(c) = datain(atlas.parcellation(c));
        plotmask(c) = roimask(atlas.parcellation(c));
    elseif isfield(atlas,'parcels')
        plotdata(c) = datain(atlas.parcels(c));
        plotmask(c) = roimask(atlas.parcels(c));
    end
    tmp = FindClosest(cindex,plotdata(c));
    if plotmask(c)
        cdata(c,:) = cmap1(tmp,:);
    else
        cdata(c,:) = cmap2(tmp,:);
    end
end

bnd.pnt = sourcemodel.pos;
bnd.tri = sourcemodel.tri;


ft_plot_mesh(bnd,'vertexcolor',cdata)
lighting gouraud
view(-90,-90)
camlight
view(0,0)

    



