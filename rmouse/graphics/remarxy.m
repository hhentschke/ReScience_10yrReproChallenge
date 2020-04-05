function remarxy(varargin)
% ** function rexy(varargin)
% resizes a 'paperposition' (=the area occupied by all axes) on a figure
% without changing its center. The idea is to minimize margins (which by
% default are a bit to spacious).
%
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% hcf               handle to axis, gca  figure to be scaled
% xfac             scalar, 1.2           paperpos will be scaled in x direction by this factor
% yfac             scalar, 1.2           paperpos will be scaled in y direction by this factor

% ----- default values & varargin -----

hcf='gcf';
% by default, make axes a trifle bigger
xfac=1.05;
yfac=1.05;

pvpmod(varargin);

if strcmp(hcf,'gcf')
  hcf=gcf;
end

oGCAUnits=get(gcf,'units');  
set(gcf,'units','normalized');  
p=get(hcf,'paperPosition');
% resize
p=[p(1)+p(3)/2*(1-xfac)  p(2)+p(4)/2*(1-yfac)  p(3)*xfac  p(4)*yfac];
set(hcf,'paperPosition',p);
set(gcf,'units',oGCAUnits);
drawnow;
