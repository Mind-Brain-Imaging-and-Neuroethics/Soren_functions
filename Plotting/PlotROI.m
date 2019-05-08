function [ROI] = PlotROI(ROI,sourcemodel,atlas)

bnd.pnt = sourcemodel.pos;
bnd.tri = sourcemodel.tri;

if ischar(ROI)
    ROI = find(strcmpi(ROI,atlas.parcellationlabel));
end
try
    plotdata = rand(8004,1)+10*(atlas.parcellation == ROI);
catch
    plotdata = rand(8004,1)+10*(atlas.parcels == ROI);
end

if isempty(which('ft_plot_mesh'))
ft_defaults
end
ft_plot_mesh(bnd,'vertexcolor',plotdata);
