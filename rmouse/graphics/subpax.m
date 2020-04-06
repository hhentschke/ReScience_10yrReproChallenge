function subpax(fh,varargin)
% ** function subpax(fh,varargin)) 
% sets x and y axes of subplots on figure to widest limits found or 
% to values in optional input variable 'xyl' (format: [x x y y]). 
% Additional input arg 'spInd' may be 
% - the index to children of type axis
% - the handles to specific axes
% Additional input arg 'ax' may be 'x', 'y' or 'xy'

xyl=[];
spInd='a';
ax='xy';
pvpmod(varargin);


if ischar(spInd) && strcmpi(spInd,'a')
  % take all
  sph=get(fh,'children');
  sph=sph(strcmp('axes',get(sph,'type')));
elseif ishandle(spInd)
  sph=spInd;
else
  % pick selected ones
  sph=get(fh,'children');
  sph=sph(strcmp('axes',get(sph,'type')));
  sph=sph(spInd);
end
  
if length(sph)>1 
  if length(xyl)==4
    xl=xyl(1:2);
    yl=xyl(3:4);
  else
    all_xl=get(sph,'xlim');
    all_yl=get(sph,'ylim');
    xl=[min([all_xl{:}]) max([all_xl{:}])];
    yl=[min([all_yl{:}]) max([all_yl{:}])];
  end
  % works vectorized!
  if ~isempty(strfind(ax,'x'))
    set(sph,'xlim',xl);
  end
  if ~isempty(strfind(ax,'y'))
    set(sph,'ylim',yl);
  end
end