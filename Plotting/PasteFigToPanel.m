function PasteFigToPanel(figfile,pasteparent,pastepanel,pastehandle,packing)
% figfile can also be a figure handle

if ischar(figfile)
   fighandle = open(figfile);
else
    fighandle = figfile;
end

figaxes = findobj('Parent',fighandle,'Type','axes');

set(0,'currentfigure',fighandle)

pastepanel.pack(packing(1),packing(2));
select
