function [ROI] = Plot_ROI_seed(ROI,corrdata,sourcemodel,atlas)

clf
subplot(2,2,1)
ROI = PlotROI(ROI,sourcemodel,atlas);

subplot(2,2,2)
corrdata = corrdata(ROI,:);

for c = 1:8004
    if isfield(atlas,'parcellation')
   plotdata(c) = corrdata(atlas.parcellation(c)); 
    else
        plotdata(c) = corrdata(atlas.parcels(c));
    end
end
if isfield(atlas,'parcellation')
plotdata(find(atlas.parcellation == ROI)) = NaN;
else
    plotdata(find(atlas.parcels == ROI)) = NaN;
end

bnd.pnt = sourcemodel.pos;
bnd.tri = sourcemodel.tri;

ft_plot_mesh(bnd,'vertexcolor',plotdata');
colorbar
clim = get(gca,'CLim');

subplot(2,2,3)

cort_size = length(find(sourcemodel.brainstructure == 1));

bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 1),:);
trindx = find(max(sourcemodel.tri,[],2) <= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);

ft_plot_mesh(bnd,'vertexcolor',plotdata(1:length(plotdata)/2)');
set(gca,'CLim',clim);


subplot(2,2,4)
bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 2),:);;
trindx = find(min(sourcemodel.tri,[],2) >= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);
bnd.tri = bnd.tri - cort_size;

ft_plot_mesh(bnd,'vertexcolor',plotdata((length(plotdata)/2+1):end)');
set(gca,'CLim',clim)
