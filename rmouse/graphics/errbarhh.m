function [axh,bh,eh]=errbarhh(x,y1,y2,e1,e2,varargin)
% ** function [axh,bh,eh]=errbarhh(x,y1,y2,e1,e2,varargin)
% produces a double horizontal bar plot where one set of bars faces towards
% the left and the other to the right; both are aligned horizontally (like
% the ones frequently used to depict the age distribution of a population).
% This is achieved by positioning two axes close to each other. Error bars
% may be plotted. Optional input arguments must be specified as parameter/
% value pairs 
%
%                         >>> INPUT VARIABLES >>>
%
% NAME              TYPE/DEFAULT          DESCRIPTION
% x                 column vector         ordinate values
% y1                array                 abscissa values for left plot
% y2                array                 abscissa values for right plot
% e1                array                 'error' values for left plot (may
%                                          be set to [])
% e2                array                 'error' values for right plot (may
%                                          be set to [])
% pos               4 element array,      position of the COMPOUND axis in
%                    [.1 .1 .8 .8]         usual matlab notation (see help
%                                          for subplot for orientation)
% barargin                                optional input arguments for barh
% invBarGroupOrd    scalar, 1             if nonzero, the order of bars
%                                          within a group will be inverted
%                                          so that the first data column in
%                                          y1 appears at the top (and not
%                                          at the bottom, which would be
%                                          the default)
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT           DESCRIPTION
% axh              2-element array        handles to axes
% bh               array                  handles to bars (barseries in
%                                          matlab V. 7.x); two columns, one
%                                          per plot
% eh               3D (!) array           handles to error bars; two slices, 
%                                          one per plot

% to do:
% - how much formatting shall be done via input args?


pos=[.1 .1 .8 .8];
barargin={};
invBarGroupOrd=0;
pvpmod(varargin);


% ----- checks of input
plotErrBar=1;
if isempty(e1) || isempty(e2)
  plotErrBar=0;
  eh1=[];
  eh2=[];
end

[n1, n2, n3]=size(y1);
if n3>1
  error('sorry, 3D not allowed');
end

if ~isequal([n1 n2], size(y2))
  error('y1 and y2 must have equal size')
end
if plotErrBar
  if ~isequal([n1 n2], size(e1)) || ~isequal([n1 n2], size(e2))
    error('e1 and e2 must have same size as y1 and y2');
  end
end
  
% ------- if requested invert bar group order
if invBarGroupOrd
  y1=fliplr(y1);
  y2=fliplr(y2);
  e1=fliplr(e1);
  e2=fliplr(e2);
end

% ------ compute positions of subplots
pos1=[pos([1 2]) pos(3)/2 pos(4)];
pos2=[pos(1)+pos(3)/2  pos(2) pos(3)/2  pos(4)];

% left bar plot: y axis on right, x axis inverted
ax1=subplot('position',pos1);
bh1=barh(x,y1,barargin{:});
set(gca,'yaxisloc','right','xdir','reverse','ytick',[]);
hold on
if plotErrBar
  % get positions of bars on axis so we can plot error bars
  for g=1:length(bh1)
    barPos = get(bh1(g),'XEndPoints')';
    eh1(g) = errorbar(y1(:, g), barPos, e1(:, g), 'horizontal', 'color', 'k', ...
        'CapSize', 0, 'linewidth', 1.5);
    eh1(g).LineStyle = 'none';
  end
  % rearrange children such that for each bar group the error lines are behind bars
  set(gca,'children',[fliplr(bh1), fliplr(eh1)]);
end
nicexyax(30);
xl=get(gca,'xlim');
set(gca,'xlim',[0 xl(2)]);
box off


% right bar plot: everything normal
ax2=subplot('position',pos2);
bh2=barh(x,y2,barargin{:});
set(gca,'ytick',[]);
hold on
if plotErrBar
  % get positions of bars on axis so we can plot error bars
  for g=1:length(bh2)
    barPos = get(bh2(g),'XEndPoints')';
    eh2(g) = errorbar(y2(:, g), barPos, e2(:, g), 'horizontal', 'color', 'k', ...
        'CapSize', 0, 'linewidth', 1.5);
    eh2(g).LineStyle = 'none';
  end
  % rearrange children such that for each bar group the error lines are behind bars
  set(gca,'children',[fliplr(bh2), fliplr(eh2)]);
end
nicexyax(30);
xl=get(gca,'xlim');
set(gca,'xlim',[0 xl(2)]);
box off

subpax(gcf);

axh=[ax1 ax2];
bh=[bh1; bh2]';
eh=cat(3,eh1,eh2);
if invBarGroupOrd
  bh=flipud(bh); 
  eh=flipdim(eh,2);
end