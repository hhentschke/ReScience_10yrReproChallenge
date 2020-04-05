function nicexax(nicefac)
if nargin<1,
  nicefac=30;
end

axis tight;
xl=get(gca,'XLim');
set(gca,'XLim',[xl(1)-(xl(2)-xl(1))/nicefac xl(2)+(xl(2)-xl(1))/nicefac]);
