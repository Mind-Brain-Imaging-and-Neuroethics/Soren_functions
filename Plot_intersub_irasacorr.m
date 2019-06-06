function [corrmats,measnames,figs] = Plot_intersub_irasacorr(outputs,measnames_bb)

if ~exist('measnames','var')

measnames_bb = {'PLE 2-8 Hz','PLE 15-85 Hz','Osci power 1.3-85 Hz','Frac power 1.3-85 Hz','OsciFrac ratio 1.3-85 Hz'};
measnames_fb = {'Delta 1.3-4 Hz','Theta 4-8 Hz','Alpha 8-13 Hz','Beta 13-30 Hz','Gamma 30-85 Hz'};
measnames_fb = [cellcat('Osci ',measnames_fb,'',0) cellcat('Frac ',measnames_fb,'',0)];

measnames = [measnames_bb measnames_fb];
end

for c = 1:length(outputs.sub)
   corrmats(:,:,c) = corr(outputs.data(:,:,c)',outputs.data(:,:,c)','Type','Spearman'); 
end

figs(1) = figure;
for c = 1:5
   subplot(2,3,c)
   imagesc(corrmats(:,:,c))
   title(measnames{c})
end
colorbar
colormap(jet)
Normalize_Clim(gcf)

figs(2) = figure;

for c = 1:10
   subplot(2,5,c)
   imagesc(corrmats(:,:,c+5))
   title(measnames{c+5})
end
colorbar
colormap(jet)
Normalize_Clim(gcf)