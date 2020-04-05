function labelscale(varargin)
% ** function labelscale(varargin)
% sets defaults for graphics such that line thickness, font size, etc. have
% a fixed default value after resizing of plots done outside matlab. An
% ordinary line plot with markers, orient portrait, A4 or the like, looks
% OK if not even a little fragile, with the default settings. ** Note: if
% figures are exported as vector graphics (e.g. postscript) and changed in
% size in e.g. Corel Draw, the thickness of the lines in the plot may not
% scale with figure size automatically - in this case ungroup elements of
% the figure and select this option manually
%
% *** copy-and-paste line:
% labelscale('scaleFac',.5,'fontSz',10,'lineW',2,'markSz',8);
% 
%                         >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT       DESCRIPTION
% scaleFac       scalar, 1.0        the factor by which the figure will be 
%                                    resized outside matlab
% fontSz         scalar, 12         the size of the font as it should 
%                                    appear AFTER resizing the plot
% lineW          scalar, 1          line width
% markSz         scalar, 8          marker size


% ----- default values & varargin -----
verbose=0;
scaleFac=1.0;
fontSz=12;
lineW=1.0;
markSz=8;

pvpmod(varargin);

set(groot,'DefaultAxesTitleFontSizeMultiplier',1.5)
% affects axes including axis labels, legends, title, etc
% set(groot,'DefaultaxesFontSize',max([4 round(fontSz/scaleFac)]));
% set(groot,'DefaultaxesLineWidth',max([.25 .25*round(2/scaleFac)]));
set(groot,'DefaultaxesFontSize',fontSz/scaleFac);
tmp1=get(groot,'DefaultaxesLineWidth');
tmp2=get(groot,'FactoryAxesTickLength');
set(groot,'DefaultAxesTickLength', tmp2+.01*tmp1);


% affects the plotted data
% set(groot,'DefaultlineLineWidth',max([.25 .25*round(4*lineW/scaleFac)]));
% set(groot,'DefaultlineMarkerSize',max([1 round(markSz/scaleFac)]));
set(groot,'DefaultlineLineWidth',lineW/scaleFac);
set(groot,'DefaultlineMarkerSize',markSz/scaleFac);

% text produced with e.g. text(...)
% set(groot,'DefaulttextFontSize',max([4 round(fontSz/scaleFac)]));
set(groot,'DefaulttextFontSize',fontSz/scaleFac);
set(groot,'DefaulttextFontWeight','normal');