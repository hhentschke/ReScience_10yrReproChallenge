function utscaleb4(xlab,ylab,varargin)
% ** function utscaleb4(xlab,ylab,varargin)
% Creates a scalebar for x and y axes of a plot of voltage versus time
% Scalebar appears in fixed position (lower left) with variable length.
%
%                         >>> INPUT VARIABLES >>>
%
% NAME              TYPE/DEFAULT          DESCRIPTION
% xlab              char arr              real unit label
% ylab              char arr, 'mV'        real unit label
% yScaleFac         array, 1              real units that correspond to a unitary 
%                                         step in the plot
% pos               char arr, 'll'        position of scale bar
%                                         'll' - lower left
%                                         'lr' - lower right
% posfac            array, [.03 .01]      fraction of x and y ranges to be added to
%                                         scalebar's position (positive values =
%                                         scale bar will be shifted outside of plot)

pos='ll';
posfac=[.03 .01];
yScaleFac=1;
pvpmod(varargin)
% this is not a good idea because the assignment of numbers to traces will be lost
%yScaleFac=unique(yScaleFac);

thickn=1.5;
xl=get(gca,'XLim');
yl=get(gca,'YLim');

% lengths (spans) of axes
axxSpan=diff(xl);
axySpan=diff(yl);

% this is the ideal fraction of the axes' length which the scale bar should have
xfrac=.1;
yfrac=.2;

% this will produce fractional lengths in [0.1..1]*[.1 .2 .4 .5 1]
xSpan=(10^(floor(log10(axxSpan))))*[.1 .2 .4 .5 1];
ySpan=(10^(floor(log10(axySpan))))*[.1 .2 .4 .5 1];
% find index to best matching fractional length
[nix,xSpanIdx]=min(abs(xSpan-xfrac*axxSpan));
[nix,ySpanIdx]=min(abs(ySpan-yfrac*axySpan));
xSpan=xSpan(xSpanIdx);
ySpan=ySpan(ySpanIdx);

% origin of scale bar on plot
switch pos
  case 'll'
    xstart=xl(1)-axxSpan*posfac(1);
    ystart=yl(1)-axySpan*posfac(2);
  case 'lr'
    xstart=xl(2)+axxSpan*posfac(1);
    ystart=yl(1)-axySpan*posfac(2);
  otherwise
    error('bad position for scalebar')
end


xtext=[sprintf('%5.3f',xSpan)  ' ' xlab];
if length(yScaleFac)>1, sstring=[repmat('%4.4g|',1,length(yScaleFac)-1) '%4.4g'];
else sstring='%4.4g';
end
ytext=[sprintf(sstring,ySpan*yScaleFac) ' ' ylab];


xtextxfac=0;
xtextyfac=0.3;
ytextxfac=0.2;
ytextyfac=0.3;

% this ensures that plot is not rescaled 
axis manual;

switch pos
  case 'll'
    lh=line([xstart xstart+xSpan xstart xstart], ...
      [ystart ystart ystart ystart+ySpan],'LineWidth', thickn, 'Color','k');
    xtextxpos=xstart;
    xtextypos=ystart-xtextyfac*ySpan;
    ytextxpos=xstart-ytextxfac*xSpan;
    ytextypos=ystart+ytextyfac*ySpan;
    rot=90;
  case 'lr'
    lh=line([xstart xstart-xSpan xstart xstart], ...
      [ystart ystart ystart ystart+ySpan],'LineWidth', thickn, 'Color','k');
    xtextxpos=xstart-xSpan;
    xtextypos=ystart-xtextyfac*ySpan;
    ytextxpos=xstart+ytextxfac*xSpan;
    ytextypos=ystart+ySpan-ytextyfac*ySpan;
    rot=-90;
end    
 
set(lh,'Clipping','off');
text(xtextxpos,xtextypos,xtext);
text(ytextxpos,ytextypos,ytext,'Rotation',rot);
