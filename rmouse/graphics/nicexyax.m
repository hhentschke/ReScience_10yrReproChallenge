function nicexyax(nicefac)
% ** function nicexyax(nicefac)
% sets x and y axis limits such that the plot 'nicely' encloses all 
% data points (input arg nicefac is optional; it is the factor determining
% relative border width; default is 20) 

if nargin<1
  nicefac=20;
end
axis tight;
xl=get(gca,'XLim');
dx=diff(xl);
yl=get(gca,'YLim');
dy=diff(yl);

if strcmp(get(gca,'xscale'),'log')
  set(gca,'XLim',[xl(1) / (1+log(xl(2)/xl(1))/nicefac), xl(2) * (1+log(xl(2)/xl(1))/nicefac)]);
else
  set(gca,'XLim',[xl(1)-dx/nicefac  xl(2)+dx/nicefac]);
end

if strcmp(get(gca,'yscale'),'log')
  set(gca,'ylim',[yl(1) / (1+log(yl(2)/yl(1))/nicefac), yl(2) * (1+log(yl(2)/yl(1))/nicefac)]);
else
  set(gca,'ylim',[yl(1)-dy/nicefac  yl(2)+dy/nicefac]);
end



