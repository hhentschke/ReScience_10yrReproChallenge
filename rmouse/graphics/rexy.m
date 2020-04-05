function rexy(varargin)
% ** function rexy(varargin)
% resizes a plot (=an axis) on a figure without changing its center.
%
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% ax               handle to axis, gca   axis to be scaled
% xfac             scalar, 1.2           axis will be scaled in x direction by this factor
% yfac             scalar, 1.2           axis will be scaled in y direction by this factor

% ----- default values & varargin -----

ax='gca';
% by default, make axes a trifle bigger
xfac=1.2;
yfac=1.2;

pvpmod(varargin);

if strcmp(ax,'gca')
  ax=gca;
end

oGCAUnits=get(gca,'units');  
set(gca,'units','normalized');  
p=get(ax,'position');
% resize
p=[p(1)+p(3)/2*(1-xfac)  p(2)+p(4)/2*(1-yfac)  p(3)*xfac  p(4)*yfac];
set(ax,'position',p);
set(gca,'units',oGCAUnits);
drawnow;
