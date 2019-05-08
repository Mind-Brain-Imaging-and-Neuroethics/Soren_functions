function PasteFig(figfile,pastehandle)
% figfile can also be a figure handle

if ischar(figfile)
   fighandle = open(figfile);
else
    fighandle = figfile;
end

figaxes = findobj('Parent',fighandle,'Type','axes');
copyobj(get(figaxes,'Children'),pastehandle)

