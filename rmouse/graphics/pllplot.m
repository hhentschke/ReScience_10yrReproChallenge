function [ylim,dy,yscaleFac,ph,xUnitFac]=pllplot(d,varargin)
%  ** function [ylim,dy,yscaleFac,ph,xUnitFac]=pllplot(d,varargin)
%   plots one or more data traces (e.g. from different channels ) on one single plot 
%   versus time, first on top, rest downwards. Each column of d represents a channel., 
%   Episodes may have different sampling intervals.
% 
% KNOWN BUGS
% - scaling (label of scale bar) is gotten wrong when 
%   -- d is a cell array of only one element 
% 
%                    >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT                DESCRIPTION
% d                2d-array or cell array      the data to be plotted
% si               scalar, 1                   sampling interval in us (all traces)
%                  array, ones(size(d,2))      sampling interval in us (for each trace individually) 
% yscale                                       determines scaling of the signals:
%                  - cell arr of strings:      data will be grouped according to strings;
%                                              within each group the channel with the max amplitude 
%                                              range will be fished out and serve as the reference 
%                                              against which all channels of this group will be normalized,
%                                              all in all resulting in a plot that is friendly to the eyes 
%                                              in terms of scaling. 
%                  - string 'si'               If yscale=='si' grouping is done according
%                                              to the sampling interval, meaning that e.g. field potential
%                                              and spike channels will be automatically up-/downscaled so
%                  - 1d-array                  amplitudes (real units) of the channels which shall have +- same amplitude on plot
% ylim             2d array                    y limits of the whole plot; if not specified, automatic setting will be applied
% noscb            scalar, 0                   if nonzero, the scalebar will be omitted (useful for multiple plots on one page)
% noplot           scalar, 0                   if nonzero, data are not plotted, but ylim, offs and yscale are returned 
%                                              (useful for multiple plots required to have the same y scales etc.)
% spacing, dy      string,scalar or 2d-array   both variables determine vertical position of traces on plot 
%                                              (dy need not be specified):
% +__________+_______________________________+_____________+________________________________________
%   SPACING  MEASURE BETWEEN                 DY             OFFSET CALCULATION
%   'maxmin' max/min of adjacent traces      []             automatic, individual offsets, spacing=1/10th of mean of amplitudes 
%                                                           of all signals
%                                            scalar         individual offsets, spacing=absolute value
%                                            2d-array
%   'fixed'  base lines (0)                  []             uniform offsets, spacing=arbitrary value
%                                            scalar         uniform offsets, spacing=scalar
%                                            2d-array       individual offsets, spacing=array
%   'percentile'                             []             0.1th and 99.9th percentile
%                                            scalar         PERCENTILE
%                                            2d-array       - not allowed -
%
% ylab           char arr, 'units'             y label for scalebar
% xOffs          scalar, 0                     x offset for plot (same unit
%                                              as si)
%
%                    <<< OUTPUT VARIABLES <<<
% NAME             TYPE/DEFAULT      DESCRIPTION
%
% ylim                               y limits of whole plot
% dy                                 spacing between individual traces on plot 
% yscaleFac         1d-array         factor by which each channel on the plot has been divided
% ph                1d-array         'plot handles' (handles to lines)
% xUnitFac          scalar           on the x axis produced by pllplot, one
%                                    unit corresponds to 1/xUnitFac
%                                    microseconds

% improvements:
% - 'order' of code, documentation
% - the documentation for spacing options is a mess and the code itself
% probably as well. 
% - allow specification of ylim as percentile of the topmost and bottommost
% channel
% - x scale should be an output variable, allowing the user to place other
% objects on the plot within proper time frame

% -------------- defaults ---------------------
% sampling interval: by default assume uniform sampling interval
noplot=0;
noscb=0;
spacing='maxmin';
yscale{1}='oS'; % oS=original scale
units=[];
verbose=0;
ylab='units';
xOffs=0;
pvpmod(varargin);
if verbose, disp(['**** ' mfilename ':']); end;

if ~isempty(units),  warndlg([mfilename ': input variable ''units'' changed to ''yscale''']); end

% -------------- determine data variable type & dimension --------------------
% cell array? 
if iscell(d)
  nchan=length(d);
  for i=1:nchan
    if ~isequal(min(size(d{i})),1),
      error('the elements of cell array d must be 1d-arrays');
    end
    ptsPChan(i)=length(d{i});
  end;
else
  % dimension of the data
  [ptsPChan,nchan]=size(d);
  % remove singleton dimensions
  d=squeeze(d);
end;
if ~all(ptsPChan)
  error('at least one channel of d is empty')
end

% -------------- initialization of variables -----------------
% check whether number of spacing values (if given for individual channels) and data channels match
if exist('dy','var'), 
  dylen=length(dy);
  if ((dylen ~= 1) & (dylen ~= nchan-1))
    warning('number of spacing values must be one or (# of channels-1) - computing ONE value for all channels');
    dy=[];
    dylen=0;
  end;  
else 
  dy=[];
  dylen=0;
end;  

% deal with sampling interval:
% - if no sampling interval was provided, use generic si=1
if ~exist('si','var')
  si=1;
  noscb=1;
end   
% - check whether number of sampling intervals (if given for individual channels) and data channels match
silen=length(si);
usi=unique(si);
if (silen ~= 1)
  if (silen ~= nchan)
    error('number of sampling interval values does not match number of channels');
  end
else  
  % - make si an array 
  si=repmat(si,1,nchan);
end;  

% ++++++++++++++  deal with units & scaling
if isstr(yscale)
  if strcmp(yscale,'si')
    yscale=cellstr([repmat('U',silen,1) num2str(makecol(si))]);
  end
end
    
% ++++++++++++++ check whether number of XXXXXXX (if given for individual channels) and data channels match
ulen=length(yscale);
if (ulen ~= 1)
  if (ulen ~= nchan)
    warning('number of unit strings does not match number of channels -  skipping scaling');
    yscale=[]; yscale{1}='oS';
  end;  
end;
[uyscale,idx1,idx2]=unique(yscale);

% yscale can be 1. cell array of strings, either one cell or as many as there are channels 2. double array
% uyscale: same, but only as many as there are unique values

% ---------------- program proper ----------------------------
% always call plotdcmprss.m first
[x,y,nsi,numOfXIntrvls,maxis,minis,xUnitLabel,xUnitFac]=plotdcmprss(d,si,xOffs,verbose);
origsi=si;
si=nsi;
rng=maxis-minis;
% in case one channel is absolutely flat rng will be 0 and may cause trouble below.
% Since this is most likely a case of a time stamp list converted into an analog stream
% assume rng=1 
rng(find(rng==0))=1;

% ---- compute/apply y scaling factors
if isnumeric(yscale)
  yscaleFac=yscale;
  uyscaleFac=uyscale;
end
% if uyscale is a one-element array or cell array there is no need to scale anything
if length(uyscale)<=1,
  if isnumeric(yscale)  
    yscaleFac=repmat(uyscale,1,nchan);    
    uyscaleFac=uyscale;    
  else
    yscaleFac=ones(1,nchan);
    uyscaleFac=1;  
  end
else
  if iscell(uyscale)
    % for the channels belonging to one type of signal (the type being identified by 
    % the strings in variable 'yscale') find mean of amplitude ranges; normalize all 
    % to an amplitude range as close to 1 as possible
    for i=1:length(uyscale)
      ugIdx=strmatch(uyscale{i},yscale);
      uyscaleFac(i)=mean(rng(ugIdx));
    end
    % adjust values in uyscale such that 'round' numbers result:
    tmpstair=[.5:.5:10]';
    tmp=repmat(tmpstair,1,length(uyscaleFac)).*10.^repmat(floor(log10(uyscaleFac)),length(tmpstair),1);
    [tmpy,tmpi]=min(abs(tmp-repmat(uyscaleFac,length(tmpstair),1)));
    uyscaleFac=tmp(sub2ind(size(tmp),tmpi,1:length(uyscaleFac)));
    clear tmp tmpi tmpy;
    yscaleFac=uyscaleFac(idx2);
  end
  % scale everything
  if iscell(y)  
    for i=1:length(yscaleFac)
      y{i}=y{i}/yscaleFac(i);
    end
  else    
    y=y./repmat(yscaleFac,size(y,1),1);
  end
  maxis=maxis./yscaleFac;
  minis=minis./yscaleFac;
  rng=maxis-minis;
end 


switch spacing
case 'maxmin'
  % if dy is anything but a scalar or an array whose length matches (number
  % of channels -1), which was checked above, try to find a reasonable
  % value based on max/min values of data
  switch dylen
    % calculate offsets between traces
  case {0,1}
    if ~dylen
      dy=mean(rng)/10.0; 
    end
    offs=cumsum([0 minis]-[maxis 0]-dy)+maxis(1)+dy;
    offs=offs(1:end-1);
  otherwise
    offs=cumsum([0 minis]-[maxis 0]-[0 dy 0])+maxis(1)+dy(1);
    offs=offs(1:end-1);
  end;
case 'fixed'
  % calculate offsets between traces
  switch dylen
  case {0,1}
    if ~ dylen, dy=-50; end;
    if dy==0
      offs=[0 0];
    else
      dy=abs(dy)*-1;
      offs=[0:dy:dy*(nchan-1)];
    end
  otherwise
    offs=cumsum([0 dy]);
  end
case 'percentile'
  % here, dy is misused as percentile value
  if dylen==0 || dylen>1
    dy=[.1 99.9];
  else
    dylen==1;
    dy=[dy 100-dy];
  end
  % compute percentiles 
  if iscell(d)
    pr=repmat(nan,2,numel(d));
    for g=1:numel(d)
      pr(1:2,g)=prctile(d{g},dy);
    end
  else
    pr=prctile(d,dy);
  end
  % flip such that upper lim is in first row
  pr=flipud(pr);
  % linearize
  pr=pr(:);
  % get rid of irrelevant values: upper lim of first channel and lower of
  % last channel
  pr([1 end])=[];
  pr=reshape(pr,2,numel(pr)/2);
  offs=cumsum([0 -diff(pr)]);
otherwise
  error('illegal spacing method required');
end

if ~exist('ylim','var')
  yl=[minis(end)+offs(end) maxis(1)+offs(1)];
  % if min of lowermost trace is more positive than max of uppermost trace
  % ylim will not contain increasing values, so just take the absolute min
  % and max
  if any(~isfinite(yl)) || diff(yl)<=0
    yl=[min([minis maxis]) max([minis maxis])];
  end
  ylim=[yl(1)-(yl(2)-yl(1))/100  yl(2)+(yl(2)-yl(1))/100];
end;  

% *** plot proper
if ~noplot,
  if iscell(y)
    for i=1:length(y)
      ph(i)=plot(x{i},y{i}+offs(i),'k-');
      hold on;
    end  
    axis tight;
  else
    y=y+repmat(offs,[size(y,1) 1]);
    ph=plot(x,y,'k-');
    axis tight;
  end;
  set(gca,'YLim',ylim);
  if ~noscb, 
    utscaleb4(xUnitLabel,ylab,'yscaleFac',uyscaleFac);
  end;
  box off; 
  axis off;
end;
dy=diff(offs);
  
% ------------------------ local function -------------------------------

function [x,y,nsi,numOfXIntrvls,maxy,miny,xUnitLabel,xUnitFac]=plotdcmprss(d,si,xOffs,verbose)
% This function is useful for getting long time series data (=tens of
% thousands of points and more) into a plot quickly. The general problem:
% Whenever the number of points to be plotted into a figure exceeds the
% horizontal screen resolution, single pixels will be plotted repeatedly
% since there are always data points that map onto one and the same pixel.
% Especially with very many data points this is a tremendous waste of
% resources.
% A solution implemented here: 
% 1. Specify the width (horizontal extent) of the figure, in pixels for
% screen plots, in dots for printouts. This is the variable 'numOfXIntrvls'
% below.
% 2. Divide data to be plotted into an according number of groups.
% 3. Within each group, determine the minimum and maximum value and retain 
% only these.
%
%                 INPUT
% Variable name             explanation
% 
% d                         data to be transformed; each column corresponds to 
%                           a channel, 3d is not allowed
% si                        the sampling interval of the data; optional
%
%
%                 OUTPUT
% y, x                      the transformed data and the 'time' vector against which it 
%                           should be plotted, respectively. Like in input variable d, 
%                           each column of y holds data from one channel, but
%                           the line order is (n denoting the index of the group)
%                             max(t(1))
%                             min(t(1))
%                             max(t(2))
%                             min(t(2))
%                               .
%                               .
%                               .
%                             max(t(n))
%                             min(t(n))                            
%                           Therefore, time vector x looks like
%                             1
%                             1
%                             2
%                             2
%                             .
%                             .      * nsi * scaling factor
%                             .
%                             n
%                             n
%                           Thus, both variables can be used as in 'plot(x,y)' 
% ..to be continued



% Offer a choice of number of intervals the original data will be divided
% into. This number should be greater than or equal to the width of the
% plot (in pixels, for screen figures; in dots, for printouts), otherwise
% the plot will look funny. Here, this number is based on the machine's
% screen resolution multiplied by a small factor to have some detail for
% zooming. It is not fixed but will be chosen from a range such that the
% division of data points by the number of intervals will result in a
% minimal remainder (0, optimally). 
set(0,'unit','pixels')
tmp=get(0,'ScreenSize');
tmp=diff(tmp([1 3]))+1;
xIntrvls=[tmp*4:tmp*6];
if ~exist('si','var'), 
  si=1; 
end

% This is the limit for the number of points per channel above which
% pllplot will switch from plotting every single data point to calling
% plotdcmprss.m which works faster with very long data traces. Since
% plotdcmprss performs quite some computations, this value should be on the 
% order of several times the resolution of the plot on the screen/printout.
% Also, it must be larger than or equal to the smallest value of xIntrvls.
pSwitchLim=xIntrvls(end)*2;


% cell array? 
if iscell(d)
  nchan=length(d);
  for i=1:nchan
    ppc(i)=length(d{i});
  end;  
  % excerpt length in ms (or points)
  excLenInMs=ppc./1000.*si(:)';
  % excerpt length in normalized units, longest=1
  excLenNorm=excLenInMs/max(excLenInMs);
  % index to longest excerpt
  [nix,longestExcIx]=max(excLenNorm);

  % determine multiplication factor for x units:
  if max(excLenInMs)<1
    % microseconds
    xUnitFac=1;
    xUnitLabel='{\mu}s';
  elseif max(excLenInMs)<1e4
    % milliseconds
    xUnitFac=1e-3;
    xUnitLabel='ms';
  else
    % seconds
    xUnitFac=1e-6;
    xUnitLabel='s';
  end    
  
  for i=1:nchan    
    % convert current channel to 1d array
    % future improvement: find all channels with same si and ppc and do
    % some of calculations below only once
    dd=d{i};
    % for each channel, do the compression trick only if the number of
    % points to be plotted is larger than or equal to the largest value of
    % xIntrvls. Don't do it with fewer points even if the density of the
    % data should be high enough to merit compression
    if ppc(i)>=pSwitchLim
      % find number of x intervals producing the least number of remaining
      % points at end
      [m,idx]=min(rem(ppc(i),round(xIntrvls/excLenNorm(i))));
      numOfXIntrvls(i)=xIntrvls(idx);
      % the number of data points falling into one interval (e.g.
      % horizontal pixel width)..
      ptsPIntrvl=floor(ppc(i)/numOfXIntrvls(i));
      % ..determines the new sampling interval 
      nsi(i)=si(i)*ptsPIntrvl;
      
      % the number of points per channel in original data which will be
      % transformed
      nppc(i)=ptsPIntrvl*numOfXIntrvls(i);
      
      % the remaining points at end 
      remYIdx=[nppc(i)+1:ppc(i)];
      if verbose
        disp([mfilename ': dividing data into ' int2str(numOfXIntrvls(i)) ' groups of ' int2str(ptsPIntrvl) ' points each, ignoring ' int2str(length(remYIdx)) ' points at end']);
      end
      
      % now divide original data into intervals, thereby creating 3d array (here with singleton dimension)
      dd=reshape(dd(1:nppc(i),:),ptsPIntrvl,numOfXIntrvls(i),1);
      
      % determine max/min and remove singleton dimension (=convert from 3d to 2d)
      tmpmaxy=max(dd,[],1);
      tmpmaxy=shiftdim(tmpmaxy,1);
      tmpminy=min(dd,[],1);
      tmpminy=shiftdim(tmpminy,1);
      
      % merge max and min values into 3d array and clear original data
      clear dd;
      tmpy=tmpmaxy;
      tmpy(:,:,2)=tmpminy;
      
      % since due to the calculations above it is no major task to
      % calculate the absolute minimum/maximum value for each channel, do
      % this now since having these values may be handy for plots of this
      % data.
      % tmpmaxy/tmpminy are row vectors holding the absolute maxima/minima
      % for each channel
      tmpmaxy=max(tmpmaxy);
      tmpminy=min(tmpminy);
      
      % reshape dd in such a way that each column again holds data from one channel, but
      % the line order is like
      % max(t(1))
      % min(t(1))
      % max(t(2))
      % min(t(2))
      %   .
      %   .
      %   .
      % max(t(n))
      % min(t(n))
      tmpy=permute(tmpy,[3 1 2]);
      tmpy=reshape(tmpy,numOfXIntrvls(i)*2,1);
      
      % a special time axis is needed:
      % 1
      % 1
      % 2
      % 2
      % .
      % .    * nsi * xUnitFac
      % .
      % n
      % n
      tmpx=repmat((1:numOfXIntrvls(i))*(nsi(i)*xUnitFac),2,1);
      tmpx=reshape(tmpx,numOfXIntrvls(i)*2,1);
    else
      % no point in compressing 
      nsi(i)=si(i);
      tmpx=((1:ppc(i))*(nsi(i)*xUnitFac))';
      tmpy=dd;
      numOfXIntrvls(i)=nan;
      tmpmaxy=max(dd);
      tmpminy=min(dd);    
    end
    x{i}=tmpx+xOffs*xUnitFac;
    y{i}=tmpy;
    maxy(i)=tmpmaxy;
    miny(i)=tmpminy;
  end
else
  % ** data d is a regular array
  if ndims(d)>2, error([mfilename ' accepts only 1d and 2d data']); end;
  % the number of points per column in original data
  ppc=size(d,1);
  % number of columns(=channels) in original data
  nchan=size(d,2);
  % length of data in ms
  excLenInMs=ppc/1000*si;

  % determine multiplication factor for x units:
  if excLenInMs<1
    % microseconds
    xUnitFac=1;
    xUnitLabel='{\mu}s';
  elseif excLenInMs<1e4
    % milliseconds
    xUnitFac=1e-3;
    xUnitLabel='ms';
  else
    % seconds
    xUnitFac=1e-6;
    xUnitLabel='s';
  end    

  if ppc>=pSwitchLim
    [m,idx]=min(rem(ppc,xIntrvls));
    numOfXIntrvls=xIntrvls(idx);
    % the number of data points falling into one interval (e.g. horizontal pixel width)..
    ptsPIntrvl=floor(ppc/numOfXIntrvls);
    % ..determines the new sampling interval 
    nsi=si*ptsPIntrvl;
    
    % the number of points per channel in original data which will be transformed
    nppc=ptsPIntrvl*numOfXIntrvls;
    
    % the remaining points at end 
    remYIdx=[nppc+1:ppc];
    if verbose
      disp([mfilename ': dividing data into ' int2str(numOfXIntrvls) ' groups of ' int2str(ptsPIntrvl) ' points each, ignoring ' int2str(length(remYIdx)) ' points at end']);
    end    
    % now divide original data into intervals, thereby creating 3d array
    d=reshape(d(1:nppc,:),ptsPIntrvl,numOfXIntrvls,nchan);
    
    % determine max/min and remove singleton dimension (=convert from 3d to 2d)
    maxy=max(d,[],1);
    maxy=shiftdim(maxy,1);
    miny=min(d,[],1);
    miny=shiftdim(miny,1);
    
    % merge max and min values into 3d array and clear original data
    clear d;
    y=maxy;
    y(:,:,2)=miny;
    
    % since due to the calculations above it is no major task to calculate the absolute 
    % minimum/maximum value for each channel, do this now since having these values may 
    % be handy for plots of this data.
    % maxy/miny are row vectors holding the absolute maxima/minima for each channel
    maxy=max(maxy);
    miny=min(miny);
    
    % reshape d in such a way that each column again holds data from one channel, but
    % the line order is like
    % max(t(1))
    % min(t(1))
    % max(t(2))
    % min(t(2))
    %   .
    %   .
    %   .
    % max(t(n))
    % min(t(n))
    y=permute(y,[3 1 2]);
    y=reshape(y,numOfXIntrvls*2,nchan);
    
    % a special time axis is needed:
    % 1
    % 1
    % 2
    % 2
    % .
    % .       * nsi * xUnitFac
    % .
    % n
    % n
    % (all nsi are identical, therefore pick just the first)
    x=repmat((1:numOfXIntrvls),2,1)*(nsi(1)*xUnitFac)+xOffs*xUnitFac; 
    x=reshape(x,numOfXIntrvls*2,1);
  else
    % no point in compressing 
    nsi=si;
    x=(1:ppc)'*(nsi*xUnitFac)+xOffs*xUnitFac;
    y=d;
    numOfXIntrvls=nan;
    maxy=max(d);
    miny=min(d);    
  end  
end