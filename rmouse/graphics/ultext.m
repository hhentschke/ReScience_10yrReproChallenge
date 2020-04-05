function th=ultext(txt,xfac,varargin)
% ** function th=ultext(txt,xfac,varargin)

% set(gca,'Clipping','off');
x=get(gca,'XLim');
y=get(gca,'YLim');
xd=get(gca,'xdir');

if ~exist('xfac','var'), xfac=0.08; end
yfac=0.92;

if strcmpi(xd,'reverse')
  th=text(x(2)-diff(x)*xfac,y(1)+diff(y)*yfac,txt,varargin{:});
else
  th=text(x(1)+diff(x)*xfac,y(1)+diff(y)*yfac,txt,varargin{:});
end
