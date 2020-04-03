function ccp01(varargin)
% generates plots of crosscorrelation functions: 
% - thgaeCC in fig. 4
% - thCC in fig. 2 
% expects r as variable input arg
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
rv={'thCCMn'};

% watch out - possible name intersection with rmouse_chan
chans={'IN 11'};
nChans=length(chans);
% x axis limits
xl=325;

% --- prepare graphics
labelscale('fontSz',8,'scaleFac',.25,'lineW',1.5,'markSz',8); 
ornt='portrait';
figdir='d:\projects\rmouse\paper_atropine\rawFig\';
printas='-dpsc2';
% printas=[];

figName=mfilename;

% --- paths, channels
rmouse_ini;
cd([WP.rootPath '\beta3_wtko\wt0001_04708']);
a001_r1;
rmouse_APcheck;
rawCh=rmouse_chan;
% locate desired channels within results structure
for i=1:nChans
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

% --- load the shit if not given as input argument
if ~exist('r','var')
  load([AP.resPath '\' AP.resFn],'r');
end

% determine index to results obtained with specified behavior in r
rix=strmatch(behav,{r(:).segmentType});
% get rid of others
r=r(rix);

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
      tmpc=reshape(tmpc,[2*AP.ccLagPts+1,nChans]);
      % cut down
      tmpc=tmpc(cix,:);
      % don't forget to flip channels more ventral than principal one  
      tmpc(:,intra_ix>AP.LFPpcInd2)=flipud(tmpc(:,intra_ix>AP.LFPpcInd2));
      c=tmpc;
      ylim=[-1.1 1.1];      
  end
  
  for chichi=1:nChans
    figure((rvi-1)*nChans+chichi); orient(ornt); 
    ph=plot(xax,c(:,chichi),'k');
    axis tight
    set(gca,'ylim',ylim);
    xl=get(gca,'xlim');
    if strcmpi(rv,'thCCMn'),
      xl=[0 xl(2)];
      set(gca,'xlim',xl);
      rexy('ax',gca,'xfac',.6,'yfac',1);      
    end
    lh=line(xl,[0 0],'linestyle','--','linewidth',get(ph,'linewidth')*.5,'color','k');
    lh=line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',get(ph,'linewidth')*.5,'color','k');    
    set(gca,'xtick',[fliplr([0:-100:xl(1)]) [100:100:xl(2)]]);
    set(gca,'ytick',[-1:.5:1]);
    if chichi==nChans
      xlabel('Lag (ms)');
    end

    if ~isempty(printas), 
      print(printas,[figdir figName '_' rv{rvi} '_' chans{chichi}]); 
    end
    
  end

%  figure(rvi); clf, orient(ornt); 
%   for chichi=1:nChans
%     plot(xax+chichi*20,c(:,chichi)+(chichi-1)*.3);
%     niceyax
%   end
%   for chichi=1:nChans
%     % center column
%     subplot(nChans,3,(chichi-1)*3+2);
%     plot(xax,c(:,chichi));
%     set(gca,'ylim',[-.5 .5]);
%     grid on
%   end
end


