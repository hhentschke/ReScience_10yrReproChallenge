function th=urtext(txt,xfac,varargin)
% ** th=urtext(txt,xfac,varargin)

set(gca,'Clipping','off');
x=get(gca,'XLim');
y=get(gca,'YLim');
xd=get(gca,'xdir');

if ~exist('xfac','var'), xfac=0.9; end
yfac=0.92;

if strcmpi(xd,'reverse')
  th=text(x(2)-diff(x)*xfac,y(1)+diff(y)*yfac,txt,varargin{:});
  set(th,'horizontalAlignment','left');    
else
  th=text(x(1)+diff(x)*xfac,y(1)+diff(y)*yfac,txt,varargin{:});
  set(th,'horizontalAlignment','right');  
end
