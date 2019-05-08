function FillBetween(xdata,ydata1,ydata2)

if size(xdata,1) > 1
   xdata = xdata'; 
end

if size(ydata1,1) > 1
   ydata1 = ydata1'; 
end

if size(ydata2,1) > 1
   ydata2 = ydata2'; 
end

x2 = [xdata, fliplr(xdata)];
inBetween = [ydata1 fliplr(ydata2)];
fill(x2,inBetween,[0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');