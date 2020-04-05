function niceyax(nicefac)
if nargin<1,
  nicefac=20;
end
axis tight;
yl=get(gca,'YLim');
set(gca,'YLim',[yl(1)-(yl(2)-yl(1))/nicefac  yl(2)+(yl(2)-yl(1))/nicefac]);
