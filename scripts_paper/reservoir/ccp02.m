function ccp02(varargin)
% generates plot of crosscorrelation functions: 
% - thgaeCC (all channels)

% ** expects r(bi) as variable input arg, bi representing behvior of
% interest **

global DS AP

if nargin>0
  r=varargin{1};
else
  butt=questdlg('loading r inside function and discarding afterwards - sure you want this?');
  if ~strcmpi(butt,'yes')
    return
  end
end

behav={'exploring'};
rv={'thgaeCCMn'};

% x axis limits
xl=250;

% --- prepare graphics
labelscale('fontSz',10,'scaleFac',.45,'lineW',1.0,'markSz',4); 
ornt='portrait';
figdir='c:\projects\hippocampus\rmouse\paper_atropine\rawFig\';
printas='-dps2';
printas=[];

figName=mfilename;

% --- paths, channels
rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
a001_r1;
rmouse_APcheck;
rawCh=rmouse_chan;

% watch out - possible name intersection with rmouse_chan
chans=AP.rawChAnNm(AP.pcInd-6:AP.pcInd)
nChan=length(chans);

% locate desired channels within results structure
for i=1:nChan
  % indices for intra-electrode CC (AP.LFPInd is not needed..)
  intra_ix(i)=strmatch(chans{i},AP.rawChAnNm(AP.LFPInd),'exact');
  % indices for inter-electrode CC
  inter_ix(i)=strmatch(chans{i},DS.rawCh(AP.allLFPIdx,1),'exact');  
  % principal channels should be clear
end
if isempty(strfind(DS.dpath,':')), DS.dpath=[WP.rootPath DS.dpath]; end

% --- all about time: si, axes
if exist([DS.dpath '\' DS.abfFn '.abf'],'file')
  abfi=abfinfo([DS.dpath '\' DS.abfFn '.abf'],'verbose',0);
  % sampling interval
  si=abfi.si;
else
  warndlg([DS.dpath '\' DS.abfFn ' does not exist; assuming sampling interval of 1000 us'])
end
xax=discrete2cont([-AP.ccLagPts:AP.ccLagPts],si*.001,'intv',0)+1;
if ~isempty(xl)
  cix=find(abs(xax)<=xl);
  xax=xax(cix);
else
  cix=1:length(xax);
end

% --- load r if not given as input argument
if ~exist('r','var')
  load([AP.resPath '\' AP.resFn],'r');
  % determine index to results obtained with specified behavior in r
  rix=strmatch(behav,{r(:).segmentType});
  % get rid of others
  r=r(rix);
end


close all
% --- extract & plot 
for rvi=1:length(rv)
  eval(['c=r.' rv{rvi} ';']);
  switch rv{rvi}
    case {'thgaeCCMn','detheCCMn'}
      c=c(cix,intra_ix);
      ylim=[-.45 .45];
    otherwise
      tmpc=[];
      % the trick: redundant combinations are empty
      tmpc=cat(1,c{setdiff(inter_ix,AP.LFPpcInd2),AP.LFPpcInd2});
      tmpc=cat(1,tmpc,c{AP.LFPpcInd2,inter_ix});
      tmpc=reshape(tmpc,[2*AP.ccLagPts+1,nChan]);
      % cut down
      tmpc=tmpc(cix,:);
      % don't forget to flip channels more ventral than principal one  
      tmpc(:,intra_ix>AP.LFPpcInd2)=flipud(tmpc(:,intra_ix>AP.LFPpcInd2));
      c=tmpc;
      ylim=[-1.1 1.1];      
  end
  
  pks=evdeal(c,'idx',{'minmaxpeak'});
  
  xoffs=0;
  yoffs=.32;
  figure(rvi); clf; orient(ornt); hold on
  % rexy('ax',gca,'xfac',.8,'yfac',1);
  % set(gca,'xaxislocation','top');
  for g=1:nChan
    ph=plot(xax+(g-1)*xoffs,c(:,nChan-g+1)+yoffs*(g-1),'k');
    ph2=plot(pks.minPeakT(nChan-g+1)+(g-1)*xoffs-xl-1,...
      pks.minPeak(nChan-g+1)+yoffs*(g-1),'ko');
    set(ph2,'markerfacecolor','k');
%     lh=line(xl,[0 0],'linestyle','--','linewidth',get(ph,'linewidth')*.5,'color','k');
    if g==nChan
      niceyax
    end
  end
  yl=get(gca,'ylim');
  xl=get(gca,'xlim');
  lh=line([0 0],yl,'linestyle','--','linewidth',get(ph,'linewidth')*.5,'color','k');
  set(gca,'xtick',[fliplr([0:-100:xl(1)]) [100:100:xl(2)]]);
  set(gca,'ytick',[1 1.5]);

%   x=repmat(xax',1,nChan);
%   y=repmat(1:nChan,length(xax),1);
%   plot3(x,y,c,'k')

  if ~isempty(printas),
    print(printas,[figdir figName '_' rv{rvi}]);
  end
end


